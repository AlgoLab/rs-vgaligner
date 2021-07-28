use ahash::RandomState;
use boomphf::hashmap::NoKeyBoomHashMap;
use bv::BitVec;
use handlegraph::hashgraph::HashGraph;
use handlegraph::handle::Edge;
use handlegraph::handlegraph::HandleGraph;
use crate::dna::reverse_complement;
use crate::kmer::{generate_kmers, generate_kmers_linearly, generate_pos_on_ref, generate_pos_on_ref_2, Kmer, KmerPos, generate_hash};
use crate::serialization::{serialize_object_to_file, deserialize_object_from_file};
use crate::utils::{find_graph_seq_length, find_forward_sequence, NodeRef};
use serde::{Deserialize, Serialize};
use rayon::prelude::ParallelSliceMut;
use std::cmp::Ordering;
use std::error::Error;

#[derive(Debug)]
pub struct Index {
    // the kmer size that this graph was built on
    pub kmer_length: u64,
    // consider only kmers where to_key(kmer) % sampling_mod == 0
    //sampling_length: u64,

    // Various data related to the starting graph -----
    // total sequence length of the graph
    seq_length: u64,
    // forward sequence of the graph, stored here for fast access during alignment
    seq_fwd: String,
    // reverse complemented sequence of the graph, for fast access during alignment
    seq_rev: String,
    // mark node starts in forward
    seq_bv: BitVec,
    // lets us map between our seq vector and handles (when our input graph is compacted!)
    //seq_by_rank: BitVec,
    // edge count
    n_edges: u64,
    // what's the most-efficient graph topology we can store?
    edges: Vec<Edge>,
    // node count
    n_nodes: u64,
    // refer to ranges in edges
    node_ref: Vec<NodeRef>,
    // End of various data related to the starting graph -----

    // Relevant part for index ---------
    // number of kmers in the index.
    // This is also the number of keys in the index
    n_kmers: u64,
    // number of kmer positions in the index.
    // Note that since a kmer can appear multiple times, n_kmer_pos >= n_kmers
    n_kmer_pos: u64,
    // our kmer hash table
    bhpf : NoKeyBoomHashMap<u64, u64>,
    // our kmer reference table (maps from bphf to index in kmer_pos_vec)
    // these are basically the hashes, aka the keys of index
    kmer_pos_ref: Vec<u64>,
    // our kmer positions table. The values of the index are positions in this vector
    kmer_pos_table: Vec<KmerPos>,
    // End of relevant part for index --------

    // if we're loaded, helps during teardown
    loaded: bool,
}

#[derive(Serialize, Deserialize)]
struct Metadata {
    seq_length : u64,
    kmer_length: u64,
    sampling_mod: f32,
    n_nodes: u64,
    n_edges: u64,
    n_kmers: u64,
    n_kmer_positions: u64
}

impl Index {

    pub fn build(
        graph: &HashGraph,
        kmer_length: u64,
        max_furcations: u64,
        max_degree: u64,
        _sampling_rate: f32,
        out_prefix: String,
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
        let mut graph_edges: Vec<Edge> = graph.edges_iter().collect();
        //let mut graph_edges: Vec<Edge> = graph.edges_par().collect();
        graph_edges.par_sort();

        // Get the forward and reverse encoded by the linearized graph
        let seq_fwd = find_forward_sequence(graph, &mut seq_bv, &mut node_ref);
        let seq_rev = reverse_complement(&seq_fwd.as_str());

        serialize_object_to_file(&seq_fwd, out_prefix.clone() + ".sqf").ok();
        serialize_object_to_file(&seq_rev, out_prefix.clone() + ".sqr").ok();
        serialize_object_to_file(&node_ref, out_prefix.clone() + ".gyn").ok();
        serialize_object_to_file(&seq_bv, out_prefix.clone() + ".sbv").ok();

        // Generate the kmers from the graph, and obtain graph-based positions.
        // Note that the resulting kmers will be already sorted according to their seq
        // and deduplicated.
        // Also, the kmers returned by generate_kmers and generate_kmers_linearly may slightly
        // differ (mostly in the forks field), but this should not impact alignment. This is due
        // to the differences in how kmers are generated by the two approaches.
        let mut kmers_on_graph: Vec<Kmer> = match graph.paths.is_empty() {
            // If paths are not available (= not provided in the input GFA)
            // use the same kmer-generation approach used in dozyg
            true => generate_kmers(
                graph,
                kmer_length as u64,
                Some(max_furcations),
                Some(max_degree),
            ),

            // Otherwise, use an optimized approach that works by exploring each path
            // linearly
            false => generate_kmers_linearly(
                graph,
                kmer_length as u64,
                Some(max_furcations),
                Some(max_degree),
            ),
        };

        serialize_object_to_file(&kmers_on_graph, out_prefix.clone() + ".kgph").ok();

        // Translate the kmers positions on the graph (obtained by the previous function)
        // into positions on the linearized forward or reverse sequence, according to which
        // strand each kmer was on.
        // This function returns:
        // - the list of kmer positions in the reference
        // - a set (actually a vec with unique values) of kmer hashes
        // - a list of offsets, that keeps track of where the positions of each kmer start
        // in the list of kmer positions
        let mut kmers_hashes : Vec<u64> = Vec::new();
        let mut kmers_start_offsets : Vec<u64> = Vec::new();
        let mut kmers_pos_on_ref : Vec<KmerPos> =
            generate_pos_on_ref_2(&graph, &kmers_on_graph, &seq_length, &node_ref,
                                  &mut kmers_hashes, &mut kmers_start_offsets);

        // TODO: add end marker between different kmers in kmers_pos_on_ref? Currently
        // it seems like looking at the next offset should be enough

        assert_eq!(kmers_hashes.len(), kmers_start_offsets.len());

        serialize_object_to_file(&kmers_hashes, out_prefix.clone() + ".kset").ok();
        serialize_object_to_file(&kmers_pos_on_ref, out_prefix.clone() + ".kpos").ok();

        // Generate a table which stores the kmers' starting offsets in a memory-efficient way,
        // as keys aren't actually stored (this is done by using a minimal perfect hash function,
        // which requires the keys to be known in advance). This however allows for false positives,
        // which will be handled in the clustering phase.
        // The table will be used to query against a new kmer in the following way:
        // - check that the new kmer is in kmers_hashes -> if not, the kmer is not in the index
        // - get the kmer's starting offset from the table (the offset refers to kmers_pos_on_ref)
        // - check all the positions in kmers_pos_on_ref until the next kmer starts
        let kmers_table =
            NoKeyBoomHashMap::new_parallel(kmers_hashes.clone(), kmers_start_offsets.clone());
        serialize_object_to_file(&kmers_table, out_prefix.clone() + ".bbx").ok();

        // Generate additional metadata, that will be used when rebuilding the index
        let meta = Metadata {
            seq_length,
            kmer_length,
            sampling_mod: _sampling_rate,
            n_nodes: graph.node_count() as u64,
            n_edges: graph.edge_count() as u64,
            n_kmers: kmers_hashes.len() as u64,
            n_kmer_positions: kmers_pos_on_ref.len() as u64
        };
        serialize_object_to_file(&meta, out_prefix.clone() + ".mtd").ok();

        // Finally, return the index
        Index {
            kmer_length,
            //sampling_length: 0,
            seq_length,
            seq_fwd,
            seq_rev,
            seq_bv,
            //seq_by_rank: Default::default(),
            n_edges: graph_edges.len() as u64,
            edges: graph_edges,
            n_nodes: number_nodes,
            node_ref,
            n_kmers: kmers_hashes.len() as u64,
            n_kmer_pos: kmers_pos_on_ref.len() as u64,
            bhpf: kmers_table,
            kmer_pos_ref: kmers_hashes,
            kmer_pos_table: kmers_pos_on_ref,
            loaded: false,
        }
    }

    pub fn load_from_prefix(out_prefix: String) -> Self {
        let seq_fwd : String = deserialize_object_from_file(out_prefix.clone() + ".sqf");
        let seq_rev : String = deserialize_object_from_file(out_prefix.clone() + ".sqr");
        let node_ref : Vec<NodeRef> = deserialize_object_from_file(out_prefix.clone() + ".gyn");
        let seq_bv : BitVec = deserialize_object_from_file(out_prefix.clone() + ".sbv");

        //let kmers_on_graph: Vec<Kmer> =
        //    deserialize_object_from_file(out_prefix.clone() + ".kgph");
        let kmer_positions_on_ref: Vec<KmerPos> =
            deserialize_object_from_file(out_prefix.clone() + ".kpos");
        let kmers_hashes : Vec<u64> =
            deserialize_object_from_file(out_prefix.clone() + ".kset");

        let kmers_table: NoKeyBoomHashMap<u64,u64> =
            deserialize_object_from_file(out_prefix.clone() + ".bbx");

        let meta : Metadata = deserialize_object_from_file(out_prefix.to_string() + ".mtd");

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

    // Find the starting position of a certain kmer in the index (or rather in kmer_pos_table)
    pub fn find_start_position_in_index(&self, seq: &str) -> Result<usize, &'static str> {
        match self.loaded {
            false => Err("Index not built yet!"),
            true => {

                if seq.len() != self.kmer_length as usize {
                    return Err("Wrong seq length, has different size from kmers")
                }

                let hash = generate_hash(&seq.to_string());

                match self.bhpf.get(&hash) {
                    Some(value) => {
                        Ok(*value as usize)
                    },
                    _ => {
                        Err("Kmer not in index")
                    }
                }

            }
        }

    }

    // Find the ending position of a certain kmer in the index (or rather in kmer_pos_table)
    pub fn find_end_position_in_index(&self, seq: &str) -> Result<usize, &'static str> {
        match self.loaded {
            false => Err("Index not built yet!"),
            true => {
                match self.find_start_position_in_index(seq) {
                    Ok(start_pos) => {

                        // Begin at the starting pos
                        let mut offset : usize = 0;
                        let mut kpos : &KmerPos = self.kmer_pos_table.get(start_pos+offset).unwrap();

                        // Step one at a time until the end is found
                        while (kpos.start != u64::max_value() && kpos.end != u64::max_value()) {
                            offset = offset + 1;
                            kpos = self.kmer_pos_table.get(start_pos+offset).unwrap();
                        }

                        // Discard delimiter, so that only the actual positions
                        // are returned
                        offset = offset - 1;

                        Ok(start_pos+offset)
                    }
                    Err(msg) => Err(msg)
                }
            }
        }

    }

    pub fn find_positions_for_query_kmer(&self, kmer : &str) -> Vec<KmerPos> {

        let mut kmer_positions_on_ref : Vec<KmerPos> = Vec::new();

        let starting_pos = self.find_start_position_in_index(kmer).unwrap();
        let ending_pos = self.find_end_position_in_index(kmer).unwrap();
        let mut offset : usize = 0;

        while starting_pos + offset <= ending_pos {
            let ref_pos : &KmerPos = self.kmer_pos_table.get(starting_pos + offset).unwrap();
            kmer_positions_on_ref.push(ref_pos.clone());
            offset += 1;
        }

        kmer_positions_on_ref

    }


}



#[cfg(test)]
mod test {
    
    use handlegraph::mutablehandlegraph::MutableHandleGraph;
    use handlegraph::pathgraph::PathHandleGraph;

    use crate::kmer::{generate_kmers_linearly, generate_hash};

    use super::*;
    use handlegraph::handle::Handle;
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

        let _forward = find_forward_sequence(&graph, &mut seq_bv, &mut node_ref);

        let kmers_graph = generate_kmers(&graph, 3, Some(100), Some(100));

        let mut kmers_hashes : Vec<u64> = Vec::new();
        let mut kmers_start_offsets : Vec<u64> = Vec::new();
        let kmers_ref = generate_pos_on_ref_2(&graph, &kmers_graph, &total_length, &node_ref, &mut kmers_hashes, &mut kmers_start_offsets);

        assert!(kmers_graph.len() < kmers_ref.len());
    }

    #[test]
    fn test_assert_both_functions_find_same_kmers() {
        let graph = create_simple_graph();

        let total_length = find_graph_seq_length(&graph);
        let mut seq_bv: BitVec = BitVec::new_fill(false, total_length + 1);
        let mut node_ref: Vec<NodeRef> = Vec::new();

        let _forward = find_forward_sequence(&graph, &mut seq_bv, &mut node_ref);

        let mut kmers_graph_dozyg = generate_kmers(&graph, 3, Some(100), Some(100));
        let mut kmers_graph_rust_ver = generate_kmers_linearly(&graph, 3, Some(100), Some(100));

        kmers_graph_dozyg.sort_by(|a,b| a.seq.cmp(&b.seq));
        kmers_graph_rust_ver.sort_by(|a,b| a.seq.cmp(&b.seq));

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
            find_forward_sequence(&graph, &mut seq_bv, &mut node_ref)
        );

        use bv::*;
        // The bitvector marks node starts, so since nodes in the linearization are
        // A -> CT -> GA -> GCA the bitvector should be...
        assert_eq!(
            bit_vec![true, true, false, true, false, true, false, false, true],
            seq_bv
        );

        assert_eq!(node_ref.len(), 4);

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
                edge_idx: 3,
                edges_to_node: 1
            }
        );
        assert_eq!(
            *node_ref.get(3).unwrap(),
            NodeRef {
                seq_idx: 5,
                edge_idx: 4,
                edges_to_node: 2
            }
        );
    }

    #[test]
    fn test_kmers_graph_generation() {
        let graph = create_simple_graph();

        let mut kmers_on_graph = generate_kmers(&graph, 3, Some(100), Some(100));

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
        let forward = find_forward_sequence(&graph, &mut seq_bv, &mut node_ref);

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
                edge_idx: 2,
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
        let seq : String = String::from("AACGT");
        let first_hash = generate_hash(&seq);
        let second_hash = generate_hash(&seq);
        assert_eq!(first_hash, second_hash);

        let different_hash = generate_hash(&String::from("AAT"));
        assert_ne!(first_hash, different_hash);
    }

    #[test]
    fn test_table() {
        let graph = create_simple_graph();

        let mut graph_edges: Vec<Edge> = graph.edges_iter().collect();
        //let mut graph_edges: Vec<Edge> = graph.edges_par().collect();
        graph_edges.par_sort();

        // Get the number of nodes in the graph
        let number_nodes = graph.graph.len() as u64;

        let seq_length = find_graph_seq_length(&graph);

        let total_length = find_graph_seq_length(&graph);
        let mut seq_bv: BitVec = BitVec::new_fill(false, total_length + 1);
        let mut node_ref: Vec<NodeRef> = Vec::new();
        let seq_fwd = find_forward_sequence(&graph, &mut seq_bv, &mut node_ref);
        let seq_rev = reverse_complement(&seq_fwd);

        let seq_fwd2 = seq_fwd.clone();
        let seq_rev2 = seq_rev.clone();

        let kmers_on_graph = generate_kmers(&graph, 3, Some(100), Some(100));

        let mut kmers_hashes : Vec<u64> = Vec::new();
        let mut kmers_start_offsets : Vec<u64> = Vec::new();
        let kmers_positions_on_ref: Vec<KmerPos> =
            generate_pos_on_ref_2(&graph, &kmers_on_graph, &seq_length, &node_ref,
                                  &mut kmers_hashes, &mut kmers_start_offsets);

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
            let ending_pos = test_index.find_end_position_in_index(&kmer.seq).unwrap();
            let mut offset : usize = 0;

            println!("Start: {}, End: {}", starting_pos, ending_pos);

            loop {
                let ref_pos : &KmerPos = kmers_positions_on_ref.get(starting_pos + offset).unwrap();

                let ref_sequence : String;
                if ref_pos.orient == true {
                    ref_sequence = seq_fwd.clone();
                } else {
                    ref_sequence = seq_rev.clone();
                }

                let ref_substring = ref_sequence.substring(ref_pos.start as usize,
                                                           ref_pos.end as usize);

                // The kmer can either be the exact substring, or a "border" of the substring
                // Since k=3, I will only make sure that at least the first and the last base are equal
                //assert_eq!(kmer.seq, ref_substring); <- not guaranteed

                println!("Kmer: {}", kmer.seq);
                println!("Ref: {}\n", ref_substring);

                assert_eq!(kmer.seq.chars().nth(0), ref_substring.chars().nth(0));
                assert_eq!(kmer.seq.chars().nth(2), ref_substring.chars().nth(ref_substring.len()-1));

                if starting_pos + offset == ending_pos {
                    break
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
        let _forward = find_forward_sequence(&graph, &mut seq_bv, &mut node_ref);

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

        let mut kmers_hashes : Vec<u64> = Vec::new();
        let mut kmers_start_offsets : Vec<u64> = Vec::new();
        let kmers_positions_on_ref: Vec<KmerPos> =
            generate_pos_on_ref_2(&graph, &kmers_on_graph, &seq_length, &node_ref,
                                  &mut kmers_hashes, &mut kmers_start_offsets);

        println!("Kmers hashes : {:#?}", kmers_hashes);
        println!("Kmers start offsets : {:#?}", kmers_start_offsets);
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
