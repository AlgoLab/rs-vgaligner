use crate::dna::reverse_complement;
use crate::kmer::{generate_hash, generate_kmers, generate_pos_on_ref_2, Kmer, KmerPos};
use crate::serialization::{deserialize_object_from_file, serialize_object_to_file};
use crate::utils::{find_forward_sequence, find_graph_seq_length, NodeRef};

use boomphf::hashmap::NoKeyBoomHashMap;
use bv::BitVec;
use handlegraph::handle::{Edge, Handle};
use handlegraph::handlegraph::HandleGraph;
use handlegraph::hashgraph::HashGraph;
use rayon::prelude::ParallelSliceMut;
use serde::{Deserialize, Serialize};

/// Represents an index over the k-mers (with k chosen as input) contained in a Handlegraph
#[derive(Debug)]
pub struct Index {
    /// The kmer size that this graph was built on
    pub kmer_length: u64,
    // consider only kmers where to_key(kmer) % sampling_mod == 0
    //sampling_length: u64,

    // Various data related to the starting graph -----
    /// Total sequence length of the graph (= length of the forward/reverse)
    seq_length: u64,
    /// Forward sequence of the graph, stored here for fast access during alignment
    pub(crate) seq_fwd: String,
    /// Reverse complemented sequence of the graph, for fast access during alignment
    pub(crate) seq_rev: String,
    /// Mark node starts in forward
    pub seq_bv: BitVec,
    // lets us map between our seq vector and handles (when our input graph is compacted!)
    //seq_by_rank: BitVec,
    /// Edge count of the input graph
    n_edges: u64,
    // what's the most-efficient graph topology we can store?
    /// Edges of the input graph
    edges: Vec<Handle>,
    /// Node count of the input graph
    n_nodes: u64,
    // refer to ranges in edges
    /// Compactly represent nodes and edges of the input graph
    node_ref: Vec<NodeRef>,
    // End of various data related to the starting graph -----

    // Relevant part for index ---------
    /// Number of (unique!) kmers in the index (this is also the number of keys in the index,
    /// also stored in [`kmer_pos_ref`])
    n_kmers: u64,
    /// Number of kmer positions in the index (note that since a kmer can appear multiple times,
    /// so the following will always hold: n_kmer_pos >= n_kmers)
    n_kmer_pos: u64,
    /// The kmer hash tablegraph_edges. This is a map that has as keys the unique kmer hashes (also found
    /// in [`kmer_pos_ref`], and as values the starting index in [`kmer_pos_table`]
    bhpf: NoKeyBoomHashMap<u64, u64>,
    // our kmer reference table (maps from bphf to index in kmer_pos_vec)
    /// The hashes over which the Index was built upon (aka the keys of [`bhpf`])
    kmer_pos_ref: Vec<u64>,
    // our kmer positions table. The values of the index are positions in this vector
    /// The oriented positions (KmerPos) in the sequence graph space. Note that, since
    /// a kmer can appear in multiple places in the graph (and also on opposite strands),
    /// it will have multiple positions in [`kmer_pos_table`]. All the positions relative
    /// to the same kmer are stored close to each other, so we only have to keep track
    /// of the starting position (this is done via [`bhpf`]), and the ending position
    /// (this is done by adding "fake" KmerPos that act as delimiters).
    kmer_pos_table: Vec<KmerPos>,
    // End of relevant part for index --------

    // if we're loaded, helps during teardown
    /// Check if the index has been loaded
    pub(crate) loaded: bool,
}

/// Additional metadata used during (de)-serialization
#[derive(Serialize, Deserialize)]
struct Metadata {
    seq_length: u64,
    kmer_length: u64,
    sampling_mod: f32,
    n_nodes: u64,
    n_edges: u64,
    n_kmers: u64,
    n_kmer_positions: u64,
}

impl Index {
    /// Build an index over the kmers of length [`kmer_length`] for a given [`graph`]
    /// and store its various files in a location with prefix [`out_prefix`].
    /// It is also possible to set a limit on the [`max_furcations`] to be used during the
    /// graph visit, and also on the [`max_degree`] that a node can have.
    pub fn build(
        graph: &HashGraph,
        kmer_length: u64,
        max_furcations: u64,
        max_degree: u64,
        _sampling_rate: f32,
        out_prefix: Option<&str>,
    ) -> Self {
        // Get the number of nodes in the graph
        let number_nodes = graph.graph.len() as u64;

        // Get the length of the sequence encoded by the graph
        let seq_length = find_graph_seq_length(graph);

        // Mark node starts in forward
        let mut seq_bv: BitVec = BitVec::new_fill(false, seq_length + 1);

        // Store offsets in fwd and edge vector
        let mut node_ref: Vec<NodeRef> = Vec::new();

        // Store edges
        let mut graph_edges: Vec<Handle> = Vec::new();
        let mut n_edges: u64 = 0;

        // Get the forward and reverse encoded by the linearized graph
        let seq_fwd = find_forward_sequence(
            graph,
            &mut seq_bv,
            &mut node_ref,
            &mut graph_edges,
            &mut n_edges,
        );
        let seq_rev = reverse_complement(&seq_fwd.as_str());

        // Generate the kmers from the graph, and obtain graph-based positions.
        // Note that the resulting kmers will be already sorted according to their seq
        // and deduplicated.
        // Also, the kmers returned by generate_kmers and generate_kmers_linearly may slightly
        // differ (mostly in the forks field), but this should not impact alignment. This is due
        // to the differences in how kmers are generated by the two approaches.
        let kmers_on_graph: Vec<Kmer> = match graph.paths.is_empty() {
            // If paths are not available (= not provided in the input GFA)
            // use the same kmer-generation approach used in dozyg
            true => generate_kmers(
                graph,
                kmer_length as u64,
                Some(max_furcations),
                Some(max_degree),
            ),

            // Otherwise, use an optimized approach that works by exploring each path
            // linearly -- currently disabled as it's bugged
            // TODO: fix this later, something is probably off with kmer equality
            // false => generate_kmers_linearly
            false => generate_kmers(
                graph,
                kmer_length as u64,
                Some(max_furcations),
                Some(max_degree),
            ),
        };

        // Translate the kmers positions on the graph (obtained by the previous function)
        // into positions on the linearized forward or reverse sequence, according to which
        // strand each kmer was on.
        // This function returns:
        // - the list of kmer positions in the reference
        // - a set (actually a vec with unique values) of kmer hashes
        // - a list of offsets, that keeps track of where the positions of each kmer start
        // in the list of kmer positions
        let mut kmers_hashes: Vec<u64> = Vec::new();
        let mut kmers_start_offsets: Vec<u64> = Vec::new();
        let kmers_pos_on_ref: Vec<KmerPos> = generate_pos_on_ref_2(
            &graph,
            &kmers_on_graph,
            &seq_length,
            &node_ref,
            &mut kmers_hashes,
            &mut kmers_start_offsets,
        );

        assert_eq!(kmers_hashes.len(), kmers_start_offsets.len());

        /*
        println!("Kmer hashes length: {}", kmers_hashes.len());
        println!("Kmer hashes: {:#?}", kmers_hashes);
        let mut new_hashes : HashSet<u64> = HashSet::from_iter(kmers_hashes.iter().cloned());
        assert_eq!(kmers_hashes.len(), new_hashes.len());
         */

        // Generate a table which stores the kmers' starting offsets in a memory-efficient way,
        // as keys aren't actually stored (this is done by using a minimal perfect hash function,
        // which requires the keys to be known in advance). This however allows for false positives,
        // which will be handled in the clustering phase.
        // The table will be used to query against a new kmer in the following way:
        // - check that the new kmer is in kmers_hashes -> if not, the kmer is not in the index
        // - get the kmer's starting offset from the table (the offset refers to kmers_pos_on_ref)
        // - check all the positions in kmers_pos_on_ref until the next kmer starts
        let kmers_table = NoKeyBoomHashMap::new_parallel(kmers_hashes.clone(), kmers_start_offsets);

        // Obtain the index
        let index = Index {
            kmer_length,
            //sampling_length: 0,
            seq_length,
            seq_fwd,
            seq_rev,
            seq_bv,
            //seq_by_rank: Default::default(),
            n_edges,
            edges: graph_edges,
            n_nodes: number_nodes,
            node_ref,
            n_kmers: kmers_hashes.len() as u64,
            n_kmer_pos: kmers_pos_on_ref.len() as u64,
            bhpf: kmers_table,
            kmer_pos_ref: kmers_hashes,
            kmer_pos_table: kmers_pos_on_ref,
            loaded: false,
        };

        println!("Index built correctly!");

        // Store the index as multiple files
        if let Some(out_prefix) = out_prefix {
            // Generate additional metadata, that will be used when rebuilding the index
            let meta = Metadata {
                seq_length,
                kmer_length,
                sampling_mod: _sampling_rate,
                n_nodes: graph.node_count() as u64,
                n_edges: graph.edge_count() as u64,
                n_kmers: index.kmer_pos_ref.len() as u64,
                n_kmer_positions: index.kmer_pos_table.len() as u64,
            };

            match index.store_with_prefix(meta, out_prefix.to_string()) {
                Err(e) => panic!("{}", e),
                _ => println!("Index stored correctly!"),
            }
        }

        index
    }

    /// Store the index in a location with prefix [out_prefix]
    fn store_with_prefix(&self, meta: Metadata, out_prefix: String) -> std::io::Result<()> {
        serialize_object_to_file(&self.seq_fwd, out_prefix.clone() + ".sqf")?;
        serialize_object_to_file(&self.seq_rev, out_prefix.clone() + ".sqr")?;
        serialize_object_to_file(&self.node_ref, out_prefix.clone() + ".gyn")?;
        serialize_object_to_file(&self.seq_bv, out_prefix.clone() + ".sbv")?;

        //serialize_object_to_file(&kmers_on_graph, out_prefix.clone() + ".kgph").ok();

        serialize_object_to_file(&self.kmer_pos_ref, out_prefix.clone() + ".kset")?;
        serialize_object_to_file(&self.kmer_pos_table, out_prefix.clone() + ".kpos")?;

        serialize_object_to_file(&self.bhpf, out_prefix.clone() + ".bbx")?;

        serialize_object_to_file(&meta, out_prefix.clone() + ".mtd")?;

        Ok(())
    }

    pub fn load_from_prefix(out_prefix: String) -> Self {
        let seq_fwd: String = deserialize_object_from_file(out_prefix.clone() + ".sqf");
        let seq_rev: String = deserialize_object_from_file(out_prefix.clone() + ".sqr");
        let node_ref: Vec<NodeRef> = deserialize_object_from_file(out_prefix.clone() + ".gyn");
        let seq_bv: BitVec = deserialize_object_from_file(out_prefix.clone() + ".sbv");

        //let kmers_on_graph: Vec<Kmer> =
        //    deserialize_object_from_file(out_prefix.clone() + ".kgph");
        let kmer_positions_on_ref: Vec<KmerPos> =
            deserialize_object_from_file(out_prefix.clone() + ".kpos");
        let kmers_hashes: Vec<u64> = deserialize_object_from_file(out_prefix.clone() + ".kset");

        let kmers_table: NoKeyBoomHashMap<u64, u64> =
            deserialize_object_from_file(out_prefix.clone() + ".bbx");

        let meta: Metadata = deserialize_object_from_file(out_prefix.to_string() + ".mtd");

        Index {
            kmer_length: meta.kmer_length,
            seq_length: meta.seq_length,
            seq_fwd,
            seq_rev,
            seq_bv,
            //seq_by_rank: Default::default(),
            n_edges: meta.n_edges,
            edges: Vec::new(),
            n_nodes: meta.n_nodes,
            node_ref,
            n_kmers: meta.n_kmers,
            n_kmer_pos: meta.n_kmer_positions,
            kmer_pos_ref: kmers_hashes,
            bhpf: kmers_table,
            kmer_pos_table: kmer_positions_on_ref,
            loaded: true,
        }
    }

    /// Find the starting position of a certain kmer with seq [`seq`] in the index
    /// (or rather in kmer_pos_table)
    fn find_start_position_in_index(&self, seq: &str) -> Result<usize, &'static str> {
        if seq.len() != self.kmer_length as usize {
            return Err("Wrong seq length, has different size from kmers");
        }

        let hash = generate_hash(&seq.to_string());

        match self.bhpf.get(&hash) {
            Some(value) => Ok(*value as usize),
            _ => Err("Kmer not in index"),
        }
    }

    /// Find the ending position of a certain kmer in the index (or rather in kmer_pos_table)
    fn find_end_position_in_index(
        &self,
        _seq: &str,
        start_pos: usize,
    ) -> Result<usize, &'static str> {
        // Begin at the starting pos
        let mut offset: usize = 0;
        let mut kpos: &KmerPos = self.kmer_pos_table.get(start_pos + offset).unwrap();

        // Step one at a time until the end is found
        while kpos.start != u64::MAX && kpos.end != u64::MAX {
            offset += 1;
            kpos = self.kmer_pos_table.get(start_pos + offset).unwrap();
        }

        // Discard delimiter, so that only the actual positions
        // are returned
        offset -= 1;

        Ok(start_pos + offset)
    }

    /// Find the positions on the graph sequence vector for the kmer
    /// having seq [`kmer_seq`].
    // NOTE: this is the equivalent of the "query" option in the original dozyg
    pub fn find_positions_for_query_kmer(&self, kmer_seq: &str) -> Vec<KmerPos> {
        let mut kmer_positions_on_ref: Vec<KmerPos> = Vec::new();

        // This also checks if the provided kmers has the same size as the
        // kmer_size used when building the Index. If that's not the case, return
        // an empty Vec
        let starting_pos = match self.find_start_position_in_index(kmer_seq) {
            Ok(pos) => pos,
            Err(_) => usize::MAX,
        };

        // Not the cleanest approach but will do for now... TODO?
        // Checking if kmer is present first would require more time...

        // Kmer is actually in the index
        if starting_pos != usize::MAX {
            let ending_pos = self
                .find_end_position_in_index(kmer_seq, starting_pos)
                .unwrap();
            let mut offset: usize = 0;

            while starting_pos + offset <= ending_pos {
                let ref_pos: &KmerPos = self.kmer_pos_table.get(starting_pos + offset).unwrap();
                kmer_positions_on_ref.push(ref_pos.clone());
                offset += 1;
            }
        }

        kmer_positions_on_ref
    }

    // ------- Bitvec operations -------

    // Get in which node a certain position is
    pub fn get_node_from_pos(&self, pos: usize, orient: bool) -> u64 {
        let rank = self.get_bv_rank(pos) as u64;

        let result = match orient {
            true => rank,
            false => self.seq_bv.len() - rank,
        };

        result
    }

    pub fn get_bv_rank(&self, pos: usize) -> usize {
        assert!(pos < self.seq_bv.len() as usize);

        let mut rank: usize = 0;
        for i in 0..pos {
            if self.seq_bv.get(i as u64) == true {
                rank += 1;
            }
        }

        rank
    }

    // Get where a certain node starts in the linearization
    pub fn get_bv_select(&self, element_no: u64) -> usize {
        let mut select: usize = 0;

        for i in 0..self.seq_bv.len() {
            if self.seq_bv.get(i) == true {
                select += 1;
            }
            if select as u64 == element_no {
                break;
            }
        }

        select
    }

    // -------
}

#[cfg(test)]
mod test {

    use handlegraph::mutablehandlegraph::MutableHandleGraph;
    use handlegraph::pathgraph::PathHandleGraph;

    use crate::kmer::{generate_hash, generate_kmers_linearly};

    use super::*;

    use substring::Substring;

    /// This function creates a simple graph, used for debugging
    ///        | 2: CT \
    /// 1:  A            4: GCA
    ///        \ 3: GA |
    fn create_simple_graph() -> HashGraph {
        let mut graph: HashGraph = HashGraph::new();

        let h1 = graph.create_handle("A".as_bytes(), 1);
        let h2 = graph.create_handle("CT".as_bytes(), 2);
        let h3 = graph.create_handle("GA".as_bytes(), 3);
        let h4 = graph.create_handle("GCA".as_bytes(), 4);

        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h1, h3));
        graph.create_edge(&Edge(h2, h4));
        graph.create_edge(&Edge(h3, h4));

        let p1 = graph.create_path_handle("P1".as_bytes(), false);
        graph.append_step(&p1, h1);
        graph.append_step(&p1, h2);
        graph.append_step(&p1, h4);

        let p2 = graph.create_path_handle("P2".as_bytes(), false);
        graph.append_step(&p2, h1);
        graph.append_step(&p2, h3);
        graph.append_step(&p2, h4);

        graph
    }

    /// This function creates a simple graph, used for debugging
    ///        | 2: T \
    ///  1: GAT         4: CA
    ///        \ 3: A |
    fn create_simple_graph_2() -> HashGraph {
        let mut graph: HashGraph = HashGraph::new();
        let h1 = graph.create_handle("GAT".as_bytes(), 1);
        let h2 = graph.create_handle("T".as_bytes(), 2);
        let h3 = graph.create_handle("A".as_bytes(), 3);
        let h4 = graph.create_handle("CA".as_bytes(), 4);

        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h1, h3));
        graph.create_edge(&Edge(h2, h4));
        graph.create_edge(&Edge(h3, h4));

        graph
    }

    #[test]
    fn test_simple_graph_2() {
        let graph = create_simple_graph_2();

        let total_length = find_graph_seq_length(&graph);
        let mut seq_bv: BitVec = BitVec::new_fill(false, total_length + 1);
        let mut node_ref: Vec<NodeRef> = Vec::new();

        let _forward =
            find_forward_sequence(&graph, &mut seq_bv, &mut node_ref, &mut vec![], &mut 0);

        let kmers_graph = generate_kmers(&graph, 3, Some(100), Some(100));

        let mut kmers_hashes: Vec<u64> = Vec::new();
        let mut kmers_start_offsets: Vec<u64> = Vec::new();
        let kmers_ref = generate_pos_on_ref_2(
            &graph,
            &kmers_graph,
            &total_length,
            &node_ref,
            &mut kmers_hashes,
            &mut kmers_start_offsets,
        );

        assert!(kmers_graph.len() < kmers_ref.len());
    }

    #[test]
    fn test_assert_both_functions_find_same_kmers() {
        let graph = create_simple_graph();

        let total_length = find_graph_seq_length(&graph);
        let mut seq_bv: BitVec = BitVec::new_fill(false, total_length + 1);
        let mut node_ref: Vec<NodeRef> = Vec::new();

        let _forward =
            find_forward_sequence(&graph, &mut seq_bv, &mut node_ref, &mut vec![], &mut 0);

        let mut kmers_graph_dozyg = generate_kmers(&graph, 3, Some(100), Some(100));
        let mut kmers_graph_rust_ver = generate_kmers_linearly(&graph, 3, Some(100), Some(100));

        kmers_graph_dozyg.sort_by(|a, b| a.seq.cmp(&b.seq));
        kmers_graph_rust_ver.sort_by(|a, b| a.seq.cmp(&b.seq));

        kmers_graph_dozyg.dedup();
        kmers_graph_rust_ver.dedup();

        //assert_eq!(kmers_graph_dozyg, kmers_graph_rust_ver);
        assert_eq!(kmers_graph_dozyg.len(), kmers_graph_rust_ver.len());

        //for kmer in &kmers_graph_rust_ver {
        //    assert!(kmers_graph_dozyg.contains(kmer));
        //}
    }

    #[test]
    fn test_forward_creation() {
        let graph = create_simple_graph();

        let total_length = find_graph_seq_length(&graph);
        let mut seq_bv: BitVec = BitVec::new_fill(false, total_length + 1);
        let mut node_ref: Vec<NodeRef> = Vec::new();

        assert_eq!(total_length, 8);
        assert_eq!(
            "ACTGAGCA",
            find_forward_sequence(&graph, &mut seq_bv, &mut node_ref, &mut vec![], &mut 0)
        );

        use bv::*;
        // The bitvector marks node starts, so since nodes in the linearization are
        // A -> CT -> GA -> GCA the bitvector should be...
        assert_eq!(
            bit_vec![true, true, false, true, false, true, false, false, true],
            seq_bv
        );

        assert_eq!(node_ref.len(), 5);

        assert_eq!(
            *node_ref.get(0).unwrap(),
            NodeRef {
                seq_idx: 0,
                edge_idx: 0,
                edges_to_node: 0
            }
        );
        assert_eq!(
            *node_ref.get(1).unwrap(),
            NodeRef {
                seq_idx: 1,
                edge_idx: 2,
                edges_to_node: 1
            }
        );
        assert_eq!(
            *node_ref.get(2).unwrap(),
            NodeRef {
                seq_idx: 3,
                edge_idx: 4,
                edges_to_node: 1
            }
        );
        assert_eq!(
            *node_ref.get(3).unwrap(),
            NodeRef {
                seq_idx: 5,
                edge_idx: 6,
                edges_to_node: 2
            }
        );
        assert_eq!(
            *node_ref.get(4).unwrap(),
            NodeRef {
                seq_idx: 8,
                edge_idx: 8,
                edges_to_node: 0
            }
        );
    }

    #[test]
    fn test_kmers_graph_generation() {
        let graph = create_simple_graph();

        let kmers_on_graph = generate_kmers(&graph, 3, Some(100), Some(100));

        assert_eq!(kmers_on_graph.len(), 14); // 8 kmers * 2 strands - 2 duplicates

        let kmers_on_graph_100 = generate_kmers(&graph, 6, Some(100), Some(100));
        assert_eq!(kmers_on_graph_100.len(), 4);

        // Check if it crashes with k=100...
        let kmers_on_graph_100 = generate_kmers(&graph, 100, Some(100), Some(100));
        assert_eq!(kmers_on_graph_100.len(), 0);
        // ...it doesn't!
    }

    #[test]
    fn test_simple_path() {
        let mut graph = HashGraph::new();

        // ACG -> TTT -> CA
        let h1 = graph.create_handle("ACG".as_bytes(), 1);
        let h2 = graph.create_handle("TTT".as_bytes(), 2);
        let h3 = graph.create_handle("CA".as_bytes(), 3);

        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h2, h3));

        let p1 = graph.create_path_handle("P1".as_bytes(), false);
        graph.append_step(&p1, h1);
        graph.append_step(&p1, h2);
        graph.append_step(&p1, h3);

        // Generate the forward
        let total_length = find_graph_seq_length(&graph);
        let mut seq_bv: BitVec = BitVec::new_fill(false, total_length + 1);
        let mut node_ref: Vec<NodeRef> = Vec::new();
        let forward =
            find_forward_sequence(&graph, &mut seq_bv, &mut node_ref, &mut vec![], &mut 0);

        assert_eq!(total_length, 8);
        assert_eq!("ACGTTTCA", forward);
        assert_eq!(
            *node_ref.get(1).unwrap(),
            NodeRef {
                seq_idx: 3,
                edge_idx: 1,
                edges_to_node: 1
            }
        );
        assert_eq!(
            *node_ref.get(2).unwrap(),
            NodeRef {
                seq_idx: 6,
                edge_idx: 3,
                edges_to_node: 1
            }
        );

        let kmers_on_graph_rust_ver = generate_kmers_linearly(&graph, 3, Some(100), Some(100));
        let kmers_on_graph_dozyg = generate_kmers(&graph, 3, Some(100), Some(100));

        assert_eq!(kmers_on_graph_rust_ver.len(), 12);
        assert_eq!(kmers_on_graph_dozyg.len(), 12);
    }

    /*
    #[test]
    fn test_self_loop() {
        let mut graph = HashGraph::new();

        //          v--
        // ACG -> TTT | -> CA

        let h1 = graph.create_handle("ACG".as_bytes(), 1);
        let h2 = graph.create_handle("TTT".as_bytes(), 2);
        let h3 = graph.create_handle("CA".as_bytes(), 3);

        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h2, h2));
        graph.create_edge(&Edge(h2, h3));

        let p1 = graph.create_path_handle("P1".as_bytes(), false);
        graph.append_step(&p1, h1);
        graph.append_step(&p1, h2);
        graph.append_step(&p1, h3);

        // Generate the forward
        let total_length = find_graph_seq_length(&graph);
        let mut seq_bv: BitVec = BitVec::new_fill(false, total_length + 1);
        let mut node_ref: Vec<NodeRef> = Vec::new();
        let forward = find_forward_sequence(&graph, &mut seq_bv, &mut node_ref);

        assert_eq!(total_length, 8);
        assert_eq!("ACGTTTCA", forward);
        assert_eq!(
            *node_ref.get(1).unwrap(),
            NodeRef {
                seq_idx: 3,
                edge_idx: 1,
                edges_to_node: 2
            }
        );
        // Note that edge_idx is +1 when compared to the previous test, i.e. there is one more edge
        // representing the loop
        assert_eq!(
            *node_ref.get(2).unwrap(),
            NodeRef {
                seq_idx: 6,
                edge_idx: 3,
                edges_to_node: 1
            }
        );

        //let kmers_on_graph_rust_ver = generate_kmers_linearly(&graph, 3, Some(100), Some(100));
        let kmers_on_graph_dozyg = generate_kmers(&graph, 3, Some(100), Some(100));

        //assert_eq!(kmers_on_graph_rust_ver.len(), 10);
        assert_eq!(kmers_on_graph_dozyg.len(), 10);

        /*
        for kmer in &kmers_on_graph_rust_ver {
            assert!(kmers_on_graph_dozyg.contains(kmer));
        }
         */
    }
     */

    #[test]
    fn test_generate_hash() {
        let seq: String = String::from("AACGT");
        let first_hash = generate_hash(&seq);
        let second_hash = generate_hash(&seq);
        assert_eq!(first_hash, second_hash);

        let different_hash = generate_hash(&String::from("AAT"));
        assert_ne!(first_hash, different_hash);
    }

    #[test]
    fn test_table() {
        let graph = create_simple_graph();

        let mut graph_edges: Vec<Handle> = Vec::new();
        let mut n_edges: u64 = 0;

        // Get the number of nodes in the graph
        let number_nodes = graph.graph.len() as u64;

        let seq_length = find_graph_seq_length(&graph);

        let total_length = find_graph_seq_length(&graph);
        let mut seq_bv: BitVec = BitVec::new_fill(false, total_length + 1);
        let mut node_ref: Vec<NodeRef> = Vec::new();
        let seq_fwd = find_forward_sequence(
            &graph,
            &mut seq_bv,
            &mut node_ref,
            &mut graph_edges,
            &mut n_edges,
        );
        let seq_rev = reverse_complement(&seq_fwd);

        let seq_fwd2 = seq_fwd.clone();
        let seq_rev2 = seq_rev.clone();

        let kmers_on_graph = generate_kmers(&graph, 3, Some(100), Some(100));

        let mut kmers_hashes: Vec<u64> = Vec::new();
        let mut kmers_start_offsets: Vec<u64> = Vec::new();
        let kmers_positions_on_ref: Vec<KmerPos> = generate_pos_on_ref_2(
            &graph,
            &kmers_on_graph,
            &seq_length,
            &node_ref,
            &mut kmers_hashes,
            &mut kmers_start_offsets,
        );

        let kmers_mphf =
            NoKeyBoomHashMap::new_parallel(kmers_hashes.clone(), kmers_start_offsets.clone());

        let test_index = Index {
            kmer_length: 3,
            //sampling_length: 0,
            seq_length,
            seq_fwd: seq_fwd2,
            seq_rev: seq_rev2,
            seq_bv,
            //seq_by_rank: Default::default(),
            n_edges: graph_edges.len() as u64,
            edges: graph_edges,
            n_nodes: number_nodes,
            node_ref,
            n_kmers: kmers_hashes.len() as u64,
            n_kmer_pos: kmers_positions_on_ref.len() as u64,
            bhpf: kmers_mphf,
            kmer_pos_ref: kmers_hashes,
            kmer_pos_table: kmers_positions_on_ref.clone(),
            loaded: true,
        };

        // Check that the hash -> start_offset mapping makes sense
        for kmer in &kmers_on_graph {
            let starting_pos = test_index.find_start_position_in_index(&kmer.seq).unwrap();
            let ending_pos = test_index
                .find_end_position_in_index(&kmer.seq, starting_pos)
                .unwrap();
            let mut offset: usize = 0;

            println!("Start: {}, End: {}", starting_pos, ending_pos);

            loop {
                let ref_pos: &KmerPos = kmers_positions_on_ref.get(starting_pos + offset).unwrap();

                let ref_sequence: String;
                if ref_pos.start_orient == true {
                    ref_sequence = seq_fwd.clone();
                } else {
                    ref_sequence = seq_rev.clone();
                }

                let ref_substring =
                    ref_sequence.substring(ref_pos.start as usize, ref_pos.end as usize);

                // The kmer can either be the exact substring, or a "border" of the substring
                // Since k=3, I will only make sure that at least the first and the last base are equal
                //assert_eq!(kmer.seq, ref_substring); <- not guaranteed

                //println!("Kmer: {}", kmer.seq);
                //println!("Ref: {}\n", ref_substring);

                assert_eq!(kmer.seq.chars().nth(0), ref_substring.chars().nth(0));
                assert_eq!(
                    kmer.seq.chars().nth(2),
                    ref_substring.chars().nth(ref_substring.len() - 1)
                );

                if starting_pos + offset == ending_pos {
                    break;
                } else {
                    offset = offset + 1;
                }
            }
        }
    }

    #[test]
    fn test_serialization() {
        let graph = create_simple_graph();

        let seq_length = find_graph_seq_length(&graph);

        // Generate the forward
        let total_length = find_graph_seq_length(&graph);
        let mut seq_bv: BitVec = BitVec::new_fill(false, total_length + 1);
        let mut node_ref: Vec<NodeRef> = Vec::new();
        let _forward =
            find_forward_sequence(&graph, &mut seq_bv, &mut node_ref, &mut vec![], &mut 0);

        // Check node_ref serialization
        let encoded_node_ref = bincode::serialize(&node_ref).unwrap();
        let deserialized_node_ref: Vec<NodeRef> = bincode::deserialize(&encoded_node_ref).unwrap();
        assert_eq!(node_ref, deserialized_node_ref);

        // Check bitvec serialization
        let encoded_seq_bv = bincode::serialize(&seq_bv).unwrap();
        let deserialized_seq_bv: BitVec = bincode::deserialize(&encoded_seq_bv).unwrap();
        assert_eq!(seq_bv, deserialized_seq_bv);

        let kmers_on_graph = generate_kmers(&graph, 3, Some(100), Some(100));

        // Check kmer_graph serialization
        let encoded_kmer_graph = bincode::serialize(&kmers_on_graph).unwrap();
        let deserialized_kmer_graph: Vec<Kmer> = bincode::deserialize(&encoded_kmer_graph).unwrap();
        assert_eq!(kmers_on_graph, deserialized_kmer_graph);

        let mut kmers_hashes: Vec<u64> = Vec::new();
        let mut kmers_start_offsets: Vec<u64> = Vec::new();
        let _kmers_positions_on_ref: Vec<KmerPos> = generate_pos_on_ref_2(
            &graph,
            &kmers_on_graph,
            &seq_length,
            &node_ref,
            &mut kmers_hashes,
            &mut kmers_start_offsets,
        );

        //println!("Kmers hashes : {:#?}", kmers_hashes);
        //println!("Kmers start offsets : {:#?}", kmers_start_offsets);
        assert_eq!(kmers_hashes.len(), kmers_start_offsets.len());

        let table: NoKeyBoomHashMap<u64, u64> =
            NoKeyBoomHashMap::new_parallel(kmers_hashes.clone(), kmers_start_offsets.clone());

        //Check kmer_table serialization
        let encoded_table = bincode::serialize(&table).unwrap();
        let deserialized_table: NoKeyBoomHashMap<u64, u64> =
            bincode::deserialize(&encoded_table).unwrap();

        for hash in kmers_hashes {
            let table_value = table.get(&hash).unwrap();
            let deserialized_table_value = deserialized_table.get(&hash).unwrap();
            assert_eq!(table_value, deserialized_table_value);
        }
    }
}
