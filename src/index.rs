use bv::BitVec;
use handlegraph::handle::{Handle,Edge,Direction,NodeId};
use handlegraph::hashgraph::HashGraph;
use handlegraph::handlegraph::HandleGraph;
use bstr::{ByteSlice, ByteVec};
use std::fs::File;
use std::io::{Write, Read};
use crate::utils::{get_bv_rank, find_sequence, NodeRef, find_sequence_po, find_graph_seq_length};
use crate::dna::reverse_complement;
use gfa::gfa::Orientation;
use crate::kmer::{generate_kmers, generate_kmers_hash, Kmer, KmerPos, generate_pos_on_forward, generate_hash, generate_mphf, merge_kmers};
use crate::io::{print_bitvec, print_seq_to_file, print_kmers_to_file, verify_kmers, verify_kmers_2, print_kmers_to_file_split};
use ahash::RandomState;
use std::path::PathBuf;
use handlegraph::mutablehandlegraph::MutableHandleGraph;
use rayon::prelude::*;

#[derive(Default)]
pub struct Index {
    // the kmer size that this graph was built on
    kmer_length : u64,
    // consider only kmers where to_key(kmer) % sampling_mod == 0
    sampling_length : u64,
    // total sequence length of the graph
    seq_length : u64,
    // forward sequence of the graph, stored here for fast access during alignment
    seq_fwd : Vec<char>,
    // reverse complemented sequence of the graph, for fast access during alignment
    seq_rev : Vec<char>,
    // mark node starts
    seq_bv : BitVec,
    // lets us map between our seq vector and handles (when our input graph is compacted!)
    seq_by_rank : BitVec, //TODO: check this
    // edge count
    n_edges : u64,
    // what's the most-efficient graph topology we can store?
    edges : Vec<Edge>,
    // node count
    n_nodes : u64,
    // refer to ranges in edges
    nodes : Vec<Handle>,
    // number of kmers in the index
    n_kmers : u64,
    // number of kmer positions in the index
    n_kmer_pos : u64,
    // our kmer hash table
    //bhpf :  TODO: hash
    // our kmer reference table (maps from bphf to index in kmer_pos_vec)
    kmer_pos_ref : Vec<u64>,
    // our kmer positions
    kmer_pos_table : Vec<u64>, //TODO: check this
    // if we're loaded, helps during teardown
    loaded : bool
}

#[derive(Default)]
pub struct Node_ref {
    seq_idx : u64,
    edge_idx : u64,
    count_prev : u64,
}


impl Index {

    //See this pattern here: https://stackoverflow.com/a/41510505/5627359
    fn new() -> Self {
        Default::default()
    }

    pub fn build(graph : &HashGraph, kmer_length : u64, max_furcations : u64, max_degree : u64,
                 sampling_rate : f32, out_prefix : &str) -> Self {

        // Create the new index
        let mut index = Index::new();

        // Get the number of nodes in the graph
        let number_nodes = graph.graph.len();

        // Get the length of the sequence encoded by the graph
        let total_length = find_graph_seq_length(graph);

        // Mark node starts in forward
        let mut seq_bv : BitVec = BitVec::new_fill(false,total_length+1);
        // Store offsets in fwd and edge vector
        let mut node_ref : Vec<NodeRef> = Vec::new();

        let forward = find_sequence_po(graph, &mut seq_bv, &mut node_ref);
        let reverse = reverse_complement(&forward.as_str());

        println!("Forward is: {}", forward);
        //println!("Reverse is: {}", reverse);
        print!("Bitvec is: ",);
        print_bitvec(&seq_bv);
        println!("\n");
        //println!("NodeRef is: {:#?}", node_ref);

        let kmers_on_graph : Vec<Kmer> = generate_kmers(graph,kmer_length as u64, Some(max_furcations), Some(max_degree));
        //println!("kmers_on_graph: {:#?}", kmers_on_graph);

        //let kmers_on_graph_rev : Vec<Kmer> = generate_kmers_rev(graph,kmer_length as u64, Some(max_furcations), Some(max_degree));
        //println!("kmers_on_graph_rev: {:#?}", kmers_on_graph_rev);

        let kmers_path_split = format!("./output/{}kmers-split.fa", out_prefix);
        //print_kmers_to_file_split(&kmers_on_graph_fwd, &kmers_on_graph_rev, &PathBuf::from(kmers_path_split));

        //let kmers_on_graph = merge_kmers(kmers_on_graph_fwd, kmers_on_graph_rev);
        //println!("kmers_on_graph: {:#?}", kmers_on_graph);

        // Store reference in a fasta file
        let reference_path = format!("./output/{}ref.fa", out_prefix);
        print_seq_to_file(&forward, &PathBuf::from(reference_path));
        // Store kmers in a fasta file
        let kmers_path = format!("./output/{}kmers.fa", out_prefix);
        print_kmers_to_file(&kmers_on_graph, &PathBuf::from(kmers_path));

        /*
        let kmers_on_seq_fwd : Vec<KmerPos> = generate_pos_on_forward(&kmers_on_graph, &forward, &seq_bv, &node_ref);
        //println!("kmers_on_seq_fwd: {:#?}", kmers_on_seq_fwd);
        let kmers_visualization_path = format!("./output/{}kmers_on_fwd.txt", out_prefix);
        verify_kmers_2(&forward, &kmers_on_graph, &kmers_on_seq_fwd, &PathBuf::from(kmers_visualization_path));

        // TODO: seed should be fixed?
        let hash_build = RandomState::new();
        let hashes = generate_hash(&kmers_on_graph, &hash_build);
        let kmer_mphf_table = generate_mphf(&hashes);
        println!("{:#?}", kmer_mphf_table);

        // Generate kmers
        let mut sorted_handles : Vec<Handle> = graph.handles_iter().collect();
        sorted_handles.sort();

        for handle in sorted_handles {
            let node = graph.get_node(&handle.id()).unwrap();
            let handle_seq = node.sequence.to_string();


            let test = graph.sequence(handle).into_string_lossy();
            println!("Forward handle {:#?} has nodeId {} with seq {}", handle, handle.unpack_number(), test);

            let test2 = graph.sequence(handle.flip()).into_string_lossy();
            let node_rev = graph.get_node(&handle.id()).unwrap();
            let handle_seq_rev = node_rev.sequence.to_string();

            println!("Backward handle {:#?} has nodeId {} with seq {}\n", handle, handle.unpack_number(), test2);
        }
         */

        index
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use handlegraph::mutablehandlegraph::MutableHandleGraph;
    use handlegraph::hashgraph::Node;
    use substring::Substring;

    /// This function creates a simple graph, used for debugging
    ///        | 2: CT \
    /// 1:  A            4: GCA
    ///        \ 3: GA |
    fn create_simple_graph() -> HashGraph {
        let mut graph : HashGraph = HashGraph::new();

        let h1 = graph.create_handle("A".as_bytes(), 1);
        let h2 = graph.create_handle("CT".as_bytes(), 2);
        let h3 = graph.create_handle("GA".as_bytes(), 3);
        let h4 = graph.create_handle("GCA".as_bytes(), 4);

        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h1, h3));
        graph.create_edge(&Edge(h2, h4));
        graph.create_edge(&Edge(h3, h4));

        graph
    }

    #[test]
    fn test_forward_creation() {
        let graph = create_simple_graph();

        let total_length = find_graph_seq_length(&graph);
        let mut seq_bv : BitVec = BitVec::new_fill(false,total_length+1);
        let mut node_ref : Vec<NodeRef> = Vec::new();

        assert_eq!(total_length, 8);
        assert_eq!("ACTGAGCA", find_sequence_po(&graph, &mut seq_bv, &mut node_ref));

        use bv::*;
        // The bitvector marks node starts, so since nodes in the linearization are
        // A -> CT -> GA -> GCA the bitvector should be...
        assert_eq!(bit_vec![true, true, false, true, false, true, false, false, true], seq_bv);

        assert_eq!(node_ref.len(), 4);
        
        assert_eq!(*node_ref.get(0).unwrap(), NodeRef { seq_idx: 0, edge_idx: 0, edges_to_node: 0 });
        assert_eq!(*node_ref.get(1).unwrap(), NodeRef { seq_idx: 1, edge_idx: 2, edges_to_node: 1 });
        assert_eq!(*node_ref.get(2).unwrap(), NodeRef { seq_idx: 3, edge_idx: 3, edges_to_node: 1 });
        assert_eq!(*node_ref.get(3).unwrap(), NodeRef { seq_idx: 5, edge_idx: 4, edges_to_node: 2 });
    }

    #[test]
    fn test_kmers_graph_generation() {
        let graph = create_simple_graph();

        let kmers_on_graph = generate_kmers(&graph, 3, Some(100), Some(100));
        println!("{:#?}", kmers_on_graph);

        assert_eq!(kmers_on_graph.len(), 12);

        let kmers_on_graph_100 = generate_kmers(&graph, 6, Some(100), Some(100));
        assert_eq!(kmers_on_graph_100.len(), 4);

        // Check if it crashes with k=100...
        let kmers_on_graph_100 = generate_kmers(&graph, 100, Some(100), Some(100));
        assert_eq!(kmers_on_graph_100.len(), 0);
        // ...it doesn't!
    }

    #[test]
    fn test_kmers_fwd_generation() {
        let graph = create_simple_graph();

        let total_length = find_graph_seq_length(&graph);
        let mut seq_bv : BitVec = BitVec::new_fill(false,total_length+1);
        let mut node_ref : Vec<NodeRef> = Vec::new();
        let forward = find_sequence_po(&graph, &mut seq_bv, &mut node_ref);

        let kmers_on_graph: Vec<Kmer> = generate_kmers(&graph, 3, Some(100), Some(100));

        // Note: forward generation and kmers_on_graph generation are independent from each other...
        // maybe they could be parallelized?

        let kmers_on_seq_fwd : Vec<KmerPos> = generate_pos_on_forward(&kmers_on_graph, &forward, &seq_bv, &node_ref);
        assert_eq!(kmers_on_graph.len(), kmers_on_graph.len());

        // There is a 1-1 relationship between kmers_on_graph and kmers_on_seq_fwd, i.e.
        // the i-eth position in each refers to the same kmer

        for i in 0..kmers_on_graph.len() {  // kmers_on_seq_fwd.len() would work as well
            let graph_kmer : &Kmer = kmers_on_graph.get(i).unwrap();
            let fwd_kmer : &KmerPos = kmers_on_seq_fwd.get(i).unwrap();

            // Check if the kmer retrieved from the graph (graph_kmer) is equal to
            // the substring of the forward starting and ending in the positions from fwd_kmer
            assert_eq!(graph_kmer.seq, forward.substring(fwd_kmer.start as usize, fwd_kmer.end as usize))
        }
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

        // Generate the forward
        let total_length = find_graph_seq_length(&graph);
        let mut seq_bv : BitVec = BitVec::new_fill(false,total_length+1);
        let mut node_ref : Vec<NodeRef> = Vec::new();
        let forward = find_sequence_po(&graph, &mut seq_bv, &mut node_ref);

        assert_eq!(total_length, 8);
        assert_eq!("ACGTTTCA", forward);
        assert_eq!(*node_ref.get(1).unwrap(),  NodeRef { seq_idx: 3, edge_idx: 1, edges_to_node: 1 });
        assert_eq!(*node_ref.get(2).unwrap(),  NodeRef { seq_idx: 6, edge_idx: 2, edges_to_node: 1 });

        let kmers_on_graph = generate_kmers(&graph, 3, Some(100), Some(100));
        assert_eq!(kmers_on_graph.len(), 6);

        let kmers_on_seq_fwd : Vec<KmerPos> = generate_pos_on_forward(&kmers_on_graph, &forward, &seq_bv, &node_ref);
    }

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

        // Generate the forward
        let total_length = find_graph_seq_length(&graph);
        let mut seq_bv : BitVec = BitVec::new_fill(false,total_length+1);
        let mut node_ref : Vec<NodeRef> = Vec::new();
        let forward = find_sequence_po(&graph, &mut seq_bv, &mut node_ref);

        assert_eq!(total_length, 8);
        assert_eq!("ACGTTTCA", forward);
        assert_eq!(*node_ref.get(1).unwrap(),  NodeRef { seq_idx: 3, edge_idx: 1, edges_to_node: 2 });
        // Note that edge_idx is +1 when compared to the previous test, i.e. there is one more edge
        // representing the loop
        assert_eq!(*node_ref.get(2).unwrap(),  NodeRef { seq_idx: 6, edge_idx: 3, edges_to_node: 1 });

        let kmers_on_graph = generate_kmers(&graph, 3, Some(100), Some(100));
        assert_eq!(kmers_on_graph.len(), 6);

    }
}