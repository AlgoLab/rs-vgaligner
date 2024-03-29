use std::ops::Range;

use boomphf::hashmap::NoKeyBoomHashMap;
use bv::BitVec;
use handlegraph::handle::{Edge, Handle};
use handlegraph::handlegraph::HandleGraph;
use handlegraph::hashgraph::HashGraph;
//use rayon::iter::IndexedParallelIterator;
//use rayon::iter::ParallelIterator;
//use rayon::prelude::IntoParallelRefIterator;
use serde::{Deserialize, Serialize};
use substring::Substring;

use crate::dna::reverse_complement;
use crate::kmer::KMER_POS_DELIMITER;
use crate::kmer::{
    generate_hash, generate_kmers, generate_kmers_parallel, generate_pos_on_ref_2, GraphKmer,
    KmerPos, SeqOrient, SeqPos,
};
use crate::serialization::{
    deserialize_object_from_file, serialize_object_to_file, store_mappings_in_file,
};
use crate::serialization::{handle_vec_deser, handle_vec_ser};
use crate::utils::{find_forward_sequence, find_graph_seq_length, NodeRef};

use crate::io::generate_json_mappings;
use log::{info, warn};
use std::time::{Duration, Instant};

/// Represents an index over the k-mers (with k chosen as input) contained in a Handlegraph.
/// This index is essentially comprised of:
/// - forward and reverse linearizations of the graph
/// - kmers represented as oriented positions (in the fwd/rev linearizations)
/// - other additional data structures that compactly represent the topology of the graph
/// (seq bv, edges, n_nodes, etc.)
#[derive(Debug, Serialize, Deserialize)]
pub struct Index {
    /// The kmer size that this graph was built on
    pub kmer_length: u64,

    // Various data related to the starting graph -----
    /// Total sequence length of the graph (= length of the forward/reverse)
    seq_length: u64,
    /// Forward sequence of the graph, stored here for fast access during alignment
    pub(crate) seq_fwd: String,
    /// Reverse complemented sequence of the graph, for fast access during alignment
    pub(crate) seq_rev: String,
    /// Mark node starts in forward
    pub seq_bv: BitVec,
    /// Edge count of the input graph
    n_edges: u64,
    /// Edges of the input graph
    #[serde(serialize_with = "handle_vec_ser")]
    #[serde(deserialize_with = "handle_vec_deser")]
    pub edges: Vec<Handle>,
    /// Node count of the input graph
    n_nodes: u64,
    /// Compactly represent nodes and edges of the input graph
    pub node_ref: Vec<NodeRef>,
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
    /// The sampling mod that was used to build the Index
    sampling_rate: Option<u64>,
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
        out_prefix: Option<&str>,
        sampling_rate: Option<u64>,
        generate_mappings: bool,
        mappings_path: Option<&str>,
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

        if generate_mappings {
            let mappings = generate_json_mappings(&graph);
            let true_mappings_path = match mappings_path {
                Some(value) => value,
                _ => "mappings.json",
            };

            match store_mappings_in_file(&mappings, true_mappings_path.to_string()) {
                Err(e) => panic!("{}", e),
                _ => info!("Mappings correctly stored in {}!", true_mappings_path),
            }
        }

        // Generate the kmers from the graph, and obtain graph-based positions.
        // Note that the resulting kmers will be already sorted according to their seq
        // and deduplicated.
        let start_generate_kmers = Instant::now();
        let kmers_on_graph: Vec<GraphKmer> = generate_kmers_parallel(
            graph,
            kmer_length as u64,
            Some(max_furcations),
            Some(max_degree),
            sampling_rate,
        );
        info!(
            "Finding the kmers required: {} ms",
            start_generate_kmers.elapsed().as_millis()
        );

        /*
        // Also, the kmers returned by generate_kmers and generate_kmers_linearly may slightly
        // differ (mostly in the forks field), but this should not impact alignment. This is due
        // to the differences in how kmers are generated by the two approaches.
        let kmers_on_graph: Vec<GraphKmer> = match graph.paths.is_empty() {
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
         */

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

        let start_convert_kmers = Instant::now();
        let kmers_pos_on_ref: Vec<KmerPos> = generate_pos_on_ref_2(
            &graph,
            kmers_on_graph,
            &seq_length,
            &node_ref,
            &mut kmers_hashes,
            &mut kmers_start_offsets,
        );
        info!(
            "Converting the kmers required: {} ms",
            start_convert_kmers.elapsed().as_millis()
        );

        assert_eq!(kmers_hashes.len(), kmers_start_offsets.len());

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
            sampling_rate,
            kmer_pos_table: kmers_pos_on_ref,
            loaded: false,
        };

        info!("Index with k={} built correctly!", index.kmer_length);
        info!(
            "Found {} different kmers, which appear in {} positions!",
            index.n_kmers, index.n_kmer_pos
        );

        // Store the index in the .idx file
        if let Some(out_prefix) = out_prefix {
            match out_prefix.ends_with(".idx") {
                false => match index.store_with_prefix(out_prefix.to_string()) {
                    Err(e) => panic!("{}", e),
                    _ => info!("Index correctly stored in {}.idx!", out_prefix.to_string()),
                },
                true => match index.store_in_file(out_prefix.to_string()) {
                    Err(e) => panic!("{}", e),
                    _ => info!("Index correctly stored in {}!", out_prefix.to_string()),
                },
            }
        }

        index
    }

    /// Store the index in a location with prefix [out_prefix]
    fn store_with_prefix(&self, out_prefix: String) -> std::io::Result<()> {
        serialize_object_to_file(&self, out_prefix.clone() + ".idx")?;
        Ok(())
    }

    /// Store the index in a file [out_file]
    fn store_in_file(&self, out_file: String) -> std::io::Result<()> {
        serialize_object_to_file(&self, out_file.clone())?;
        Ok(())
    }

    /// Load the index in a location with prefix [out_prefix]
    pub fn load_from_prefix(out_prefix: String) -> Self {
        let index: Index = deserialize_object_from_file(out_prefix.to_string() + ".idx");
        index
    }

    /// Load the index from a file [out_file]
    pub fn load_from_file(out_file: String) -> Self {
        let index: Index = deserialize_object_from_file(out_file.to_string());
        index
    }

    /// Find the starting position of a certain kmer with seq [`seq`] in the index
    /// (or rather in kmer_pos_table)
    fn find_start_position_in_index(&self, seq: &str) -> Result<usize, &'static str> {
        if seq.len() != self.kmer_length as usize {
            return Err("Wrong seq length, has different size from kmers");
        }

        let hash = generate_hash(&seq.to_string());

        match self.sampling_rate.is_none()
            || (self.sampling_rate.is_some() && hash % self.sampling_rate.unwrap() == 0)
        {
            true => match self.kmer_pos_ref.contains(&hash) {
                true => Ok(*self.bhpf.get(&hash).unwrap() as usize),
                false => Err("Kmer not in index"),
            },
            false => Err("Hash not in sampling rate"),
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
        while *kpos != KMER_POS_DELIMITER {
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

    /// Find to which node (of the original graph) a certain position (in the linearization)
    /// belongs to.
    pub fn node_id_from_seqpos(&self, pos: &SeqPos) -> u64 {
        //let rank = self.get_bv_rank(pos.position as usize) as u64;
        //println!("fwd: {:#?}", self.seq_fwd);
        //println!("rev: {:#?}", self.seq_rev);
        //println!("Bitvec is: {:#?}, rank: {}", self.seq_bv, rank);
        //assert!(rank < self.node_ref.len() as u64 - 1);

        // IMPORTANT NOTE: no matter the orientation, the returned
        // node id should always be the same. The orient only
        // matters when actually creating the Handle, where
        // even => forward, odd => reverse.
        let result = match pos.orient {
            SeqOrient::Forward => {
                self.get_bv_rank(pos.position as usize) as u64
            },
            SeqOrient::Reverse => {
                // TODO: to +1 or not to +1?
                self.n_nodes - self.get_bv_inverse_rank(pos.position as usize) as u64 + 1
            } //self.n_nodes - rank + 1,
            //SeqOrient::Reverse => rank,
        };

        result
    }

    /// Find to which handle (of the original graph) a certain position (in the linearization)
    /// belongs to.
    pub fn handle_from_seqpos(&self, pos: &SeqPos) -> Handle {
        let node_id = self.node_id_from_seqpos(pos);
        //println!("Node is: {}", node_id);
        let handle = match pos.orient {
            SeqOrient::Forward => Handle::from_integer(node_id * 2),
            SeqOrient::Reverse => Handle::from_integer(node_id * 2 + 1),
        };
        handle
    }

    /// Find the rank (starting from the left) of a given position.
    /// Useful for finding the nodeid of a seqpos with the forward orient.
    pub fn get_bv_rank(&self, pos: usize) -> usize {
        assert!(pos < self.seq_bv.len() as usize);

        let mut rank: usize = 0;

        for i in 0..=pos {
            if self.seq_bv.get(i as u64) == true {
                rank += 1;
            }
        }

        rank
    }

    /// Find the rank (starting from the left) of a given position.
    /// Useful for finding the nodeid of a seqpos with the reverse orient.
    pub fn get_bv_inverse_rank(&self, pos: usize) -> usize {
        assert!(pos < self.seq_bv.len() as usize);

        let mut inverse_rank: usize = 0;

        // -1 because last 1 (=end marker) does not matter
        let start_point = (self.seq_bv.len() - 1) as usize;

        for i in (start_point - pos..=start_point).rev() {
            if self.seq_bv.get(i as u64) == true {
                inverse_rank += 1;
            }
        }

        inverse_rank
    }

    // Get where a certain node starts in the linearization
    pub fn get_bv_select(&self, element_no: u64) -> u64 {
        if element_no == 0 {
            panic!("Element_no should be > 0");
        }

        let mut select: usize = 0;
        let mut start_pos: u64 = 0;

        for i in 0..self.seq_bv.len() {
            if self.seq_bv.get(i) == true {
                select += 1;
            }
            if select as u64 == element_no {
                start_pos = i;
                break;
            }
        }

        start_pos
    }

    // -------

    // ------- Various access methods for our Index -------

    // TODO: maybe add some checks for input handle (not really necessary but...)

    /// Find the position of a certain handle in the noderef.
    pub fn noderef_pos_from_handle(&self, handle: &Handle) -> usize {
        (u64::from(handle.id()) - 1) as usize
    }

    /// Find the noderef of a certain handle.
    pub fn noderef_from_handle(&self, handle: &Handle) -> &NodeRef {
        //Obtain id and remove 1 (because kmerpos are 0-based)
        let node_ref_pos = self.noderef_pos_from_handle(handle);
        self.node_ref.get(node_ref_pos).unwrap()
    }

    /// Find the sequence (=label) of a given handle.
    /// This is the equivalent of graph.sequence(handle), but it does
    /// not require the original graph.
    pub fn seq_from_handle(&self, handle: &Handle) -> String {
        let curr_node_ref = self.noderef_from_handle(handle);
        let curr_handle_pos = self.noderef_pos_from_handle(handle);

        // Last node is a marker
        assert!(
            curr_handle_pos < self.node_ref.len() - 1,
            "Handle is: {:#?}, nodeid: {}",
            handle,
            handle.id()
        );

        let next_node_ref = self.node_ref.get(curr_handle_pos + 1).unwrap();
        let ref_seq: String = match handle.is_reverse() {
            false => self.seq_fwd.clone(),
            _ => self.seq_rev.clone(),
        };

        let (start, end): (usize, usize) = match handle.is_reverse() {
            false => (
                curr_node_ref.seq_idx as usize,
                next_node_ref.seq_idx as usize,
            ),
            _ => (
                ref_seq.len() - next_node_ref.seq_idx as usize,
                ref_seq.len() - curr_node_ref.seq_idx as usize,
            ),
        };

        ref_seq.substring(start, end).to_string()
    }

    /// Find all the edges that start/end in a given handle.
    pub fn edges_from_handle(&self, handle: &Handle) -> &[Handle] {
        assert!(u64::from(handle.id()) <= self.n_nodes);

        let edges_interval = self.edges_interval_from_handle(handle);
        let edges: &[Handle] = &self.edges[edges_interval];
        edges
    }

    /// Find the interval (in index.edges) of the edges for a given handle.
    fn edges_interval_from_handle(&self, handle: &Handle) -> Range<usize> {
        assert!(u64::from(handle.id()) <= self.n_nodes);

        let node_ref_pos = (u64::from(handle.id()) - 1) as usize;

        let node_ref = self.node_ref.get(node_ref_pos).unwrap();
        let next_node_ref = self.node_ref.get(node_ref_pos + 1).unwrap();

        node_ref.edge_idx as usize..next_node_ref.edge_idx as usize
    }

    /// Find the incoming edges to a given handle.
    /// This is the equivalent of graph.edges_iter() with a left orient, but
    /// does not require the graph.
    pub fn incoming_edges_from_handle(&self, handle: &Handle) -> Vec<Handle> {
        assert!(u64::from(handle.id()) <= self.n_nodes);

        let edges_interval = self.edges_interval_from_handle(handle);
        let node_ref = self.noderef_from_handle(handle);

        let incoming_edges: Vec<Handle> = match handle.is_reverse() {
            false => {
                let incoming_edges_interval =
                    edges_interval.start..edges_interval.start + node_ref.edges_to_node as usize;
                self.edges[incoming_edges_interval].to_vec()
            }
            true => self
                .outgoing_edges_from_handle(&handle.flip())
                .iter()
                .map(|x| x.flip())
                .rev()
                .collect(),
        };
        incoming_edges
    }

    /// Find the outgoing edges of a given handle.
    /// This is the equivalent of graph.edges_iter() with a right orient, but
    /// does not require the graph.
    pub fn outgoing_edges_from_handle(&self, handle: &Handle) -> Vec<Handle> {
        assert!(u64::from(handle.id()) <= self.n_nodes);

        let edges_interval = self.edges_interval_from_handle(handle);
        let node_ref = self.noderef_from_handle(handle);

        let outgoing_edges: Vec<Handle> = match handle.is_reverse() {
            false => {
                let outgoing_edges_interval =
                    edges_interval.start + node_ref.edges_to_node as usize..edges_interval.end;
                self.edges[outgoing_edges_interval].to_vec()
            }
            true => self
                .incoming_edges_from_handle(&handle.flip())
                .iter()
                .map(|x| x.flip())
                .rev()
                .collect::<Vec<Handle>>()
                .to_vec(),
        };

        outgoing_edges
    }

    /// Find the sequence implied by a range of seqpos.
    pub fn seq_from_start_end_seqpos(&self, begin: &SeqPos, end: &SeqPos) -> String {
        let substring = match (begin.orient, end.orient) {
            (SeqOrient::Forward, SeqOrient::Forward) => self
                .seq_fwd
                .substring(begin.position as usize, end.position as usize)
                .to_string(),
            (SeqOrient::Reverse, SeqOrient::Reverse) => self
                .seq_rev
                .substring(begin.position as usize, end.position as usize)
                .to_string(),
            // TODO: this is 100% not right, maybe I should take a part from fwd and another from rev?
            _ => self
                .seq_fwd
                .substring(begin.position as usize, end.position as usize)
                .to_string(),
        };
        substring
    }
    // -------
}

#[cfg(test)]
mod test {
    use bstr::ByteVec;
    use handlegraph::mutablehandlegraph::MutableHandleGraph;
    use handlegraph::pathgraph::PathHandleGraph;
    use itertools::Itertools;
    use substring::Substring;

    use crate::kmer::{generate_hash, generate_kmers_linearly, get_seq_pos, SeqOrient, SeqPos};

    use super::*;
    use crate::io::QuerySequence;
    use gfa::gfa::GFA;
    use gfa::parser::GFAParser;
    use std::path::PathBuf;

    /// This function creates a simple graph, used for debugging
    ///          | 2: CT \
    /// FWD  1: A         4: GCA
    ///          \ 3: GA |
    ///
    ///          | 2: AG \
    /// REV  1: T         4: TGC
    ///          \ 3: TC |
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
    ///             | 2: T \
    ///  FWD:  1: GAT         4: CA
    ///             \ 3: A |
    ///
    ///             | 2: A \
    ///  REV:  1: ATC        4: TG
    ///             \ 3: T |
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
            kmers_graph,
            &total_length,
            &node_ref,
            &mut kmers_hashes,
            &mut kmers_start_offsets,
        );

        //assert!(kmers_graph.len() < kmers_ref.len());
    }

    /*
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
     */

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

        // Only for testing
        let kmers_on_graph_clone = kmers_on_graph.clone();

        let mut kmers_hashes: Vec<u64> = Vec::new();
        let mut kmers_start_offsets: Vec<u64> = Vec::new();
        let kmers_positions_on_ref: Vec<KmerPos> = generate_pos_on_ref_2(
            &graph,
            kmers_on_graph,
            &seq_length,
            &node_ref,
            &mut kmers_hashes,
            &mut kmers_start_offsets,
        );
        println!("Kmer pos on ref: {:#?}", kmers_positions_on_ref);

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
            sampling_rate: None,
            kmer_pos_table: kmers_positions_on_ref.clone(),
            loaded: true,
        };

        // Check that the hash -> start_offset mapping makes sense
        for kmer in &kmers_on_graph_clone {
            let starting_pos = test_index.find_start_position_in_index(&kmer.seq).unwrap();
            let ending_pos = test_index
                .find_end_position_in_index(&kmer.seq, starting_pos)
                .unwrap();
            let mut offset: usize = 0;

            println!("Start: {}, End: {}", starting_pos, ending_pos);

            loop {
                let ref_pos: &KmerPos = kmers_positions_on_ref.get(starting_pos + offset).unwrap();

                let ref_sequence: String;
                if ref_pos.start.orient == SeqOrient::Forward {
                    ref_sequence = seq_fwd.clone();
                } else {
                    ref_sequence = seq_rev.clone();
                }

                let ref_substring = ref_sequence.substring(
                    ref_pos.start.position as usize,
                    ref_pos.end.position as usize,
                );

                // The kmer can either be the exact substring, or a "border" of the substring
                // Since k=3, I will only make sure that at least the first and the last base are equal
                //assert_eq!(kmer.seq, ref_substring); <- not guaranteed

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
        let index = Index::build(&graph, 3, 100, 100, None, None, false, None);
        let serialized_index = bincode::serialize(&index).unwrap();
        let deserialized_index: Index = bincode::deserialize(&serialized_index).unwrap();

        assert_eq!(index.kmer_length, deserialized_index.kmer_length);
        assert_eq!(index.seq_length, deserialized_index.seq_length);
        assert_eq!(index.seq_fwd, deserialized_index.seq_fwd);
        assert_eq!(index.seq_rev, deserialized_index.seq_rev);
        assert_eq!(index.seq_bv, deserialized_index.seq_bv);
        assert_eq!(index.n_edges, deserialized_index.n_edges);
        assert_eq!(index.edges, deserialized_index.edges);
        assert_eq!(index.n_nodes, deserialized_index.n_nodes);
        assert_eq!(index.node_ref, deserialized_index.node_ref);
        assert_eq!(index.n_kmers, deserialized_index.n_kmers);
        assert_eq!(index.n_kmer_pos, deserialized_index.n_kmer_pos);
        assert_eq!(index.kmer_pos_ref, deserialized_index.kmer_pos_ref);
        assert_eq!(index.kmer_pos_table, deserialized_index.kmer_pos_table);
        assert_eq!(index.loaded, deserialized_index.loaded);

        /*
        for hash in kmers_hashes {
            let table_value = table.get(&hash).unwrap();
            let deserialized_table_value = deserialized_table.get(&hash).unwrap();
            assert_eq!(table_value, deserialized_table_value);
        }
         */
    }

    #[test]
    fn test_index_access() {
        // Build the index
        let graph = create_simple_graph();
        let index = Index::build(&graph, 3, 100, 100, None, None, false, None);

        // find the positions for "ACT" kmer
        let ref_pos_vec = index.find_positions_for_query_kmer("ACT");
        assert_eq!(ref_pos_vec.len(), 1);

        let act_pos = KmerPos {
            start: SeqPos {
                orient: SeqOrient::Forward,
                position: 0,
            },
            end: SeqPos {
                orient: SeqOrient::Forward,
                position: 3,
            },
        };
        assert_eq!(*ref_pos_vec.get(0).unwrap(), act_pos)
    }

    #[test]
    fn test_index_access_2() {
        // Build the index
        let mut graph = HashGraph::new();

        // TTT -> AAA
        let h1 = graph.create_handle("TTT".as_bytes(), 1);
        let h2 = graph.create_handle("AAA".as_bytes(), 2);
        graph.create_edge(&Edge(h1, h2));

        let index = Index::build(&graph, 3, 100, 100, None, None, false, None);

        // find the positions for "ACT" kmer
        let ref_pos_vec = index.find_positions_for_query_kmer("TTT");
        assert_eq!(ref_pos_vec.len(), 2);

        let ttt_pos_fwd = KmerPos {
            start: SeqPos {
                orient: SeqOrient::Forward,
                position: 0,
            },
            end: SeqPos {
                orient: SeqOrient::Forward,
                position: 3,
            },
        };
        assert_eq!(*ref_pos_vec.get(0).unwrap(), ttt_pos_fwd);

        let ttt_pos_rev = KmerPos {
            start: SeqPos {
                orient: SeqOrient::Reverse,
                position: 0,
            },
            end: SeqPos {
                orient: SeqOrient::Reverse,
                position: 3,
            },
        };
        assert_eq!(*ref_pos_vec.get(1).unwrap(), ttt_pos_rev)
    }

    #[test]
    fn test_index_access_bv() {
        // Build the index
        let graph = create_simple_graph();
        let index = Index::build(&graph, 3, 100, 100, None, None, false, None);

        //println!("Node ref: {:#?}", index.node_ref);
        for handle in graph.handles_iter().sorted() {
            //println!("Handle: {:#?}, id: {}", handle, (u64::from(handle.id())-1) as usize);

            //Obtain id and remove 1 (this is necessary to access the correct node ref)
            let node_ref_pos = (u64::from(handle.id()) - 1) as usize;

            let node_ref = index.node_ref.get(node_ref_pos).unwrap();
            println!(
                "Handle: {:#?}, Nodeid: {}, Node_ref: {:#?}",
                handle,
                handle.id(),
                node_ref
            );
            println!("REV Handle: {:#?}", handle.flip());
        }
    }

    #[test]
    fn test_index_access_edges() {
        // Build the index
        let graph = create_simple_graph();
        let index = Index::build(&graph, 3, 100, 100, None, None, false, None);

        for handle in graph.handles_iter().sorted() {
            let node_ref_pos = (u64::from(handle.id()) - 1) as usize;

            let node_ref = index.node_ref.get(node_ref_pos).unwrap();
            let next_node_ref = index.node_ref.get(node_ref_pos + 1).unwrap();

            let edges_idx = node_ref.edge_idx;
            let next_edges_idx = next_node_ref.edge_idx;

            let edges: &[Handle] = &index.edges[edges_idx as usize..next_edges_idx as usize];
            for edge in edges {
                println!("Node: {}, Edge: {:#?}", handle.id(), edge);
            }
        }
    }

    #[test]
    fn test_index_access_nodes() {
        // Build the index
        let graph = create_simple_graph();
        let index = Index::build(&graph, 3, 100, 100, None, None, false, None);

        let pos1 = SeqPos {
            orient: SeqOrient::Forward,
            position: 0,
        };
        assert_eq!(index.node_id_from_seqpos(&pos1), 1); //First node_id is 1 even in the graph

        //println!("BV: {:#?}", index.seq_bv);
        //println!("FWD: {:#?}", index.seq_fwd);
        let pos2 = SeqPos {
            orient: SeqOrient::Forward,
            position: 2,
        };
        assert_eq!(index.node_id_from_seqpos(&pos2), 2);

        let pos3 = SeqPos {
            orient: SeqOrient::Reverse,
            position: 0,
        };
        assert_eq!(index.node_id_from_seqpos(&pos3), 4);
    }

    #[test]
    fn test_compare_sequential_parallel_graphkmer() {
        let graph = create_simple_graph();
        let kmers_on_graph = generate_kmers(&graph, 3, Some(100), Some(100));
        let kmers_on_graph_parallel =
            generate_kmers_parallel(&graph, 3, Some(100), Some(100), None);
        assert_eq!(kmers_on_graph, kmers_on_graph_parallel);

        let graph2 = create_simple_graph_2();
        let kmers_on_graph2 = generate_kmers(&graph2, 3, Some(100), Some(100));
        let kmers_on_graph_parallel2 =
            generate_kmers_parallel(&graph2, 3, Some(100), Some(100), None);
        assert_eq!(kmers_on_graph2, kmers_on_graph_parallel2)
    }

    #[test]
    fn test_edges_from_handle() {
        // Build the index
        let graph = create_simple_graph();
        let index = Index::build(&graph, 3, 100, 100, None, None, false, None);
        let mut handles: Vec<Handle> = graph.handles_iter().collect();
        handles.sort();

        assert_eq!(
            index.edges_from_handle(handles.get(0).unwrap()),
            &vec![*handles.get(1).unwrap(), *handles.get(2).unwrap()]
        );
        assert_eq!(
            index.edges_from_handle(handles.get(1).unwrap()),
            &vec![*handles.get(0).unwrap(), *handles.get(3).unwrap()]
        );
        assert_eq!(
            index.edges_from_handle(handles.get(2).unwrap()),
            &vec![*handles.get(0).unwrap(), *handles.get(3).unwrap()]
        );
        assert_eq!(
            index.edges_from_handle(handles.get(3).unwrap()),
            &vec![*handles.get(1).unwrap(), *handles.get(2).unwrap()]
        );
    }

    #[test]
    fn test_index_incoming_outgoing_edges() {
        // Build the index
        let graph = create_simple_graph();
        let index = Index::build(&graph, 3, 100, 100, None, None, false, None);

        let mut handles: Vec<Handle> = graph.handles_iter().collect();
        handles.sort();
        println!("Handles: {:#?}", handles);
        println!("Edges: {:#?}", index.edges);
        println!("NodeRef: {:#?}", index.node_ref);

        assert_eq!(
            index.incoming_edges_from_handle(handles.get(0).unwrap()),
            vec![]
        );
        assert_eq!(
            index.outgoing_edges_from_handle(handles.get(0).unwrap()),
            vec![*handles.get(1).unwrap(), *handles.get(2).unwrap()]
        );

        assert_eq!(
            index.incoming_edges_from_handle(handles.get(1).unwrap()),
            vec![*handles.get(0).unwrap()]
        );
        assert_eq!(
            index.outgoing_edges_from_handle(handles.get(1).unwrap()),
            vec![*handles.get(3).unwrap()]
        );

        assert_eq!(
            index.incoming_edges_from_handle(handles.get(2).unwrap()),
            vec![*handles.get(0).unwrap()]
        );
        assert_eq!(
            index.outgoing_edges_from_handle(handles.get(2).unwrap()),
            vec![*handles.get(3).unwrap()]
        );

        assert_eq!(
            index.incoming_edges_from_handle(handles.get(3).unwrap()),
            vec![*handles.get(1).unwrap(), *handles.get(2).unwrap()]
        );
        assert_eq!(
            index.outgoing_edges_from_handle(handles.get(3).unwrap()),
            vec![]
        );

        // REVERSE
        assert_eq!(
            index.incoming_edges_from_handle(&handles.get(0).unwrap().flip()),
            vec![
                handles.get(2).unwrap().flip(),
                handles.get(1).unwrap().flip()
            ]
        );
        assert_eq!(
            index.outgoing_edges_from_handle(&handles.get(0).unwrap().flip()),
            vec![]
        );

        assert_eq!(
            index.incoming_edges_from_handle(&handles.get(3).unwrap().flip()),
            vec![]
        );
        assert_eq!(
            index.outgoing_edges_from_handle(&handles.get(3).unwrap().flip()),
            vec![
                handles.get(2).unwrap().flip(),
                handles.get(1).unwrap().flip()
            ]
        );

        assert_eq!(
            index.incoming_edges_from_handle(&handles.get(1).unwrap().flip()),
            vec![handles.get(3).unwrap().flip()]
        );
        assert_eq!(
            index.outgoing_edges_from_handle(&handles.get(1).unwrap().flip()),
            vec![handles.get(0).unwrap().flip()]
        );
    }

    #[test]
    fn test_index_seq_from_start_end_seqpos() {
        // Build the index
        let graph = create_simple_graph();
        let index = Index::build(&graph, 3, 100, 100, None, None, false, None);

        let begin = SeqPos::new(SeqOrient::Forward, 0);
        let end = SeqPos::new(SeqOrient::Forward, index.seq_length);
        assert_eq!(index.seq_from_start_end_seqpos(&begin, &end), index.seq_fwd);

        let begin2 = SeqPos::new(SeqOrient::Reverse, 0);
        let end2 = SeqPos::new(SeqOrient::Reverse, index.seq_length);
        assert_eq!(
            index.seq_from_start_end_seqpos(&begin2, &end2),
            index.seq_rev
        );

        let begin3 = SeqPos::new(SeqOrient::Forward, 0);
        let end3 = SeqPos::new(SeqOrient::Forward, 3);
        assert_eq!(
            index.seq_from_start_end_seqpos(&begin3, &end3),
            "ACT".to_string()
        );
    }

    #[test]
    fn test_seq_from_handle() {
        // Build the index
        let graph = create_simple_graph();
        let index = Index::build(&graph, 3, 100, 100, None, None, false, None);
        let mut handles: Vec<Handle> = graph.handles_iter().collect();
        handles.sort();

        //println!("FWD: {:#?} REV: {:#?}", index.seq_fwd, index.seq_rev);
        for handle in handles {
            assert_eq!(
                graph.sequence(handle).into_string_lossy(),
                index.seq_from_handle(&handle)
            );
            let rev = handle.clone().flip();
            println!(
                "ID is: {} Handle is: {:#?}, ID Rev is: {} Rev handle is: {:#?}",
                handle.id(),
                handle,
                rev.id(),
                rev
            );
            //println!("FWD seq: {:#?}, REV seq: {:#?}", graph.sequence(handle).into_string_lossy(), graph.sequence(rev).into_string_lossy());
            assert_eq!(
                graph.sequence(rev).into_string_lossy(),
                index.seq_from_handle(&rev)
            )
        }
    }

    #[test]
    fn test_handle_from_seqpos() {
        // Build the index
        let graph = create_simple_graph();
        let index = Index::build(&graph, 3, 100, 100, None, None, false, None);
        let mut handles: Vec<Handle> = graph.handles_iter().sorted().collect();

        let p1: SeqPos = SeqPos::new(SeqOrient::Forward, 0);
        assert_eq!(index.handle_from_seqpos(&p1), *handles.get(0).unwrap());

        for handle in graph.handles_iter().sorted() {
            println!("Handle is {:#?}", handle);
        }

        let p2: SeqPos = SeqPos::new(SeqOrient::Reverse, 0);
        assert_eq!(
            index.handle_from_seqpos(&p2),
            handles.last().unwrap().flip()
        );
    }

    #[test]
    fn test_reverse_handles() {
        let mut graph = HashGraph::new();
        let h1 = graph.create_handle("AAA".as_bytes(), 1);
        let h2 = graph.create_handle("TTT".as_bytes(), 2);
        let h3 = graph.create_handle("CCC".as_bytes(), 3);
        let h4 = graph.create_handle("GGG".as_bytes(), 4);

        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h1, h3));
        graph.create_edge(&Edge(h2, h4));
        graph.create_edge(&Edge(h3, h4));

        let index = Index::build(&graph, 3, 100, 100, None, None, false, None);

        let mut handles: Vec<Handle> = graph.handles_iter().collect();
        handles.sort();
        for fwd_handle in handles {
            let rev_handle = fwd_handle.flip();
            let rev_seq = graph.sequence(rev_handle).into_string_lossy();
            let _fwd_seq = graph.sequence(fwd_handle).into_string_lossy();
            //println!("Id: {}({}), Fwd is: {:#?}, Rev is: {:#?}", fwd_handle.id(), rev_handle.id(), fwd_seq, rev_seq);
            let kpos = index.find_positions_for_query_kmer(rev_seq.as_str());
            for pos in kpos {
                //println!("Pos is: {:#?}", pos.start);
                let retrieved_handle = index.handle_from_seqpos(&pos.start);
                if retrieved_handle.is_reverse() {
                    assert_eq!(index.handle_from_seqpos(&pos.start), rev_handle);
                }
            }
        }
    }

    #[test]
    fn test_seqpos_returns_all() {
        let graph = create_simple_graph();
        let index = Index::build(&graph, 3, 100, 100, None, None, false, None);
        assert_eq!(index.seq_fwd.len(), index.seq_rev.len());
        for i in 0..index.seq_fwd.len() {
            for orient in vec![SeqOrient::Forward, SeqOrient::Reverse] {
                index.handle_from_seqpos(&SeqPos::new(orient, i as u64));
            }
        }
    }

    #[test]
    fn test_wrong_index() {
        let mut graph = HashGraph::new();
        let h1 = graph.create_handle("AAAAAAA".as_bytes(), 1);
        let h2 = graph.create_handle("TTT".as_bytes(), 2);
        let h3 = graph.create_handle("CCC".as_bytes(), 3);
        let h4 = graph.create_handle("GGGGGGG".as_bytes(), 4);
        let h5 = graph.create_handle("GGG".as_bytes(), 5);
        let h6 = graph.create_handle("CCC".as_bytes(), 6);
        let h7 = graph.create_handle("TTTTTTT".as_bytes(), 7);

        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h1, h3));
        graph.create_edge(&Edge(h2, h4));
        graph.create_edge(&Edge(h3, h4));
        graph.create_edge(&Edge(h4, h5));
        graph.create_edge(&Edge(h4, h6));
        graph.create_edge(&Edge(h5, h7));
        graph.create_edge(&Edge(h6, h7));

        let seq_length = find_graph_seq_length(&graph);

        let mut seq_bv: BitVec = BitVec::new_fill(false, seq_length + 1);
        let mut node_ref: Vec<NodeRef> = Vec::new();
        let mut graph_edges: Vec<Handle> = Vec::new();
        let mut n_edges: u64 = 0;
        let seq_fwd = find_forward_sequence(
            &graph,
            &mut seq_bv,
            &mut node_ref,
            &mut graph_edges,
            &mut n_edges,
        );
        let seq_rev = reverse_complement(&seq_fwd.as_str());

        let kmer_length = 11;
        let max_furcations = 100;
        let max_degree = 100;
        let kmers_on_graph: Vec<GraphKmer> = generate_kmers_parallel(
            &graph,
            kmer_length as u64,
            Some(max_furcations),
            Some(max_degree),
            None,
        );
        //println!("graph k: {:#?}", kmers_on_graph);

        let mut sorted_handles: Vec<Handle> = graph.handles_iter().sorted().collect();

        for handle in sorted_handles {
            println!("Handle is: {:#?}, id: {}", handle, u64::from(handle.id()));
            println!(
                "seqpos: {}",
                get_seq_pos(&handle, &node_ref, &seq_length, &1)
            );
            let flipped = handle.flip();
            println!("Handle flipped? {} : {:#?}", flipped.is_reverse(), flipped);
            println!(
                "seqpos: {}",
                get_seq_pos(&flipped, &node_ref, &seq_length, &0)
            );
            println!();
        }

        let mut kmers_hashes: Vec<u64> = Vec::new();
        let mut kmers_start_offsets: Vec<u64> = Vec::new();
        let kmers_pos_on_ref: Vec<KmerPos> = generate_pos_on_ref_2(
            &graph,
            kmers_on_graph,
            &seq_length,
            &node_ref,
            &mut kmers_hashes,
            &mut kmers_start_offsets,
        );

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
            n_nodes: graph.handles_iter().count() as u64,
            node_ref,
            n_kmers: kmers_hashes.len() as u64,
            n_kmer_pos: kmers_pos_on_ref.len() as u64,
            bhpf: kmers_table,
            kmer_pos_ref: kmers_hashes,
            sampling_rate: None,
            kmer_pos_table: kmers_pos_on_ref,
            loaded: false,
        };

        for handle in graph.handles_iter().sorted() {
            println!("Id is: {}", handle.unpack_number());
        }

        for i in 0..index.node_ref.len() - 1 {
            let curr_node_ref = index.node_ref.get(i).unwrap();
            assert_eq!(
                (i + 1) as u64,
                u64::from(index.node_id_from_seqpos(&SeqPos {
                    orient: SeqOrient::Forward,
                    position: curr_node_ref.seq_idx
                }))
            );
            assert_eq!(
                (index.n_nodes as usize - i) as u64,
                u64::from(index.node_id_from_seqpos(&SeqPos {
                    orient: SeqOrient::Reverse,
                    position: curr_node_ref.seq_idx
                }))
            );
        }

        /*
        let query = QuerySequence {
            name: "".to_string(),
            seq: "CACGAGAGAAACTCCCCCGTGATTCCTAATACTGGGAGTCAAAGAGAACTGCTCATCAGTTCATCTGAAGGATGGAATCTCCAGCGAGACTAGATTCCTT".to_string()
        };

        let query_kmers = query.split_into_kmers(kmer_length as usize);
        let mut positions : Vec<KmerPos> = Vec::new();
        for kmer in query_kmers {
            positions.extend(index.find_positions_for_query_kmer(kmer.as_str()).into_iter());
            //assert!(kmers_on_graph.iter().any(|k| k.seq == kmer))
        }

        for pos in positions {
            println!("Handle: {:#?}, Node id: {}", index.handle_from_seqpos(&pos.start).as_integer(), u64::from(index.handle_from_seqpos(&pos.start).id()));
        }

        assert!(positions
            .iter()
            .any(|p| u64::from(index.handle_from_seqpos(&p.start).id()) == 10));
         */
    }

    #[test]
    fn test_inverse_rank() {
        let graph = create_simple_graph();
        let index = Index::build(&graph, 3, 100, 100, None, None, false, None);

        println!("Seq bv: {:#?}", index.seq_bv);
        let ranks: Vec<usize> = (0..index.seq_bv.len() - 1)
            .into_iter()
            .map(|i| index.get_bv_rank(i as usize))
            .collect();
        let inverse_ranks: Vec<usize> = (0..index.seq_bv.len() - 1)
            .into_iter()
            .map(|i| index.get_bv_inverse_rank(i as usize))
            .collect();

        assert_eq!(ranks, vec![1, 2, 2, 3, 3, 4, 4, 4]);
        assert_eq!(inverse_ranks, vec![1, 1, 1, 2, 2, 3, 3, 4]);
    }

    #[test]
    fn test_index_returns_same_positions() {
        // Build the index
        let graph = create_simple_graph();
        let index = Index::build(&graph, 3, 100, 100, None, None, false, None);

        for handle in graph.handles_iter().sorted() {
            let node_start_bv = index.get_bv_select(u64::from(handle.id()));

            let handle_rank = handle.unpack_number() - 1;
            let node_start_noderef = index.node_ref.get(handle_rank as usize).unwrap().seq_idx;

            assert_eq!(node_start_bv as u64, node_start_noderef);
        }
    }

    #[test]
    fn test_index_contains_multinode_kmers() {
        // Build the index
        let graph = create_simple_graph();
        let index = Index::build(&graph, 5, 100, 100, None, None, false, None);

        assert!(index.find_positions_for_query_kmer(&"ACTGC").len() > 0);
        assert!(index.find_positions_for_query_kmer(&"CTGCA").len() > 0);

        let mut graph2: HashGraph = HashGraph::new();
        let h1 = graph2.create_handle("ACG".as_bytes(), 1);
        let h2 = graph2.create_handle("C".as_bytes(), 2);
        let h3 = graph2.create_handle("G".as_bytes(), 3);
        let h4 = graph2.create_handle("TTTTT".as_bytes(), 4);
        graph2.create_edge(&Edge(h1, h2));
        graph2.create_edge(&Edge(h1, h3));
        graph2.create_edge(&Edge(h2, h4));
        graph2.create_edge(&Edge(h3, h4));

        let index2 = Index::build(&graph2, 5, 100, 100, None, None, false, None);

        let acggt_pos = index2.find_positions_for_query_kmer(&"ACGGT");
        assert!(acggt_pos.len() > 0);
        let first_pos = acggt_pos.get(0).unwrap();
        assert_eq!(first_pos.start.position, 0);
        assert_eq!(first_pos.end.position, 6);

        let gcttt_pos = index2.find_positions_for_query_kmer(&"GCTTT");
        assert!(gcttt_pos.len() > 0);
        let first_pos = gcttt_pos.get(0).unwrap();
        assert_eq!(first_pos.start.position, 2);
        assert_eq!(first_pos.end.position, 8);

        let ctttt_pos = index2.find_positions_for_query_kmer(&"CTTTT");
        assert!(ctttt_pos.len() > 0);
        let first_pos = ctttt_pos.get(0).unwrap();
        assert_eq!(first_pos.start.position, 3);
        assert_eq!(first_pos.end.position, 9);

        let mut graph3: HashGraph = HashGraph::new();
        let h1 = graph3.create_handle("ACG".as_bytes(), 1);
        let h2 = graph3.create_handle("C".as_bytes(), 2);
        let h3 = graph3.create_handle("G".as_bytes(), 3);
        let h4 = graph3.create_handle("TTTTT".as_bytes(), 4);
        let h5 = graph3.create_handle("TA".as_bytes(), 5);
        let h6 = graph3.create_handle("CG".as_bytes(), 6);
        let h7 = graph3.create_handle("TTT".as_bytes(), 7);

        graph3.create_edge(&Edge(h1, h2));
        graph3.create_edge(&Edge(h1, h3));
        graph3.create_edge(&Edge(h2, h4));
        graph3.create_edge(&Edge(h3, h4));
        graph3.create_edge(&Edge(h4, h5));
        graph3.create_edge(&Edge(h4, h6));
        graph3.create_edge(&Edge(h5, h7));
        graph3.create_edge(&Edge(h6, h7));

        let index3 = Index::build(&graph3, 5, 100, 100, None, None, false, None);

        let ttcgt_pos = index3.find_positions_for_query_kmer(&"TTCGT");
        assert!(ttcgt_pos.len() > 0);
        let first_pos = ttcgt_pos.get(0).unwrap();
        assert_eq!(first_pos.start.position, 8);
        assert_eq!(first_pos.end.position, 15);
    }
}
