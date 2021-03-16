use bv::BitVec;
use handlegraph::handle::{Handle,Edge,Direction,NodeId};
use handlegraph::hashgraph::HashGraph;
use handlegraph::handlegraph::HandleGraph;
use bstr::ByteSlice;
use std::fs::File;
use std::io::{Write, Read};
use crate::utils::{get_bv_rank, find_sequence, NodeRef, find_sequence_po};
use crate::dna::reverse_complement;
use gfa::gfa::Orientation;
use crate::kmer::{generate_kmers, generate_kmers_hash};

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

// TODO: add remaining traits etc.

/*pub trait IndexUtilities {
    fn build(graph : &HashGraph,
             kmer_length : &u64,
             max_furcations : &u64,
             max_degree : &u64,
             sampling_rate : &f32,
             out_prefix : &str) -> Index;
    //fn load() -> Index;
}*/

impl Index {

    //See this pattern here: https://stackoverflow.com/a/41510505/5627359
    fn new() -> Self {
        Default::default()
    }

    pub fn build(graph : &HashGraph, kmer_length : u64, max_furcations : u64, max_degree : u64,
                 sampling_rate : f32, out_prefix : &str) -> Self {

        // Create the new index
        let mut index = Index::new();

        let number_nodes = graph.graph.len();

        // Find total length of the sequences in the graph
        let mut total_length = 0;
        for value in graph.handles_iter() {
            let node = graph.get_node(&value.id()).unwrap();
            total_length += node.sequence.len() as u64;
        }
        index.seq_length = total_length;

        // Mark node starts in forward
        let mut seq_bv : BitVec = BitVec::new_fill(false,total_length+1);
        // Store offsets in fwd and edge vector
        let mut node_ref : Vec<NodeRef> = Vec::new();

        let forward = find_sequence_po(graph, &mut seq_bv, &mut node_ref);
        let reverse = reverse_complement(&forward.as_str());

        println!("Forward is: {}", forward);
        println!("Reverse is: {}", reverse);
        //println!("BV is: {:#?}", seq_bv);
        //println!("Node_ref is: {:#?}", node_ref);

        let kmers_on_graph = generate_kmers(graph,kmer_length as u64, Some(max_degree));
        println!("kmers_on_graph: {:#?}", kmers_on_graph);

        //let kmers_on_seq_fwd : Vec<KmerSeq> = generate_pos_on_fwd(kmers_on_graph, seq_bv, node_ref);
        //let hashes = generate_kmers_hash(&kmers_on_seq_fwd);
        //println!("hashes: {:#?}", hashes);

        index
    }
}