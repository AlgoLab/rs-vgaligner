use bv::BitVec;
use handlegraph::handle::{Handle,Edge,Direction,NodeId};
use handlegraph::hashgraph::HashGraph;
use handlegraph::handlegraph::HandleGraph;
use bstr::ByteSlice;
use std::fs::File;
use std::io::{Write, Read};
use crate::utils::{get_bv_rank, find_sequence, NodeRef, find_sequence_po, find_graph_seq_length};
use crate::dna::reverse_complement;
use gfa::gfa::Orientation;
use crate::kmer::{generate_kmers, generate_kmers_hash, Kmer, KmerPos, generate_pos_on_forward, generate_hash, generate_mphf};
use crate::io::{print_bitvec, print_seq_to_file, print_kmers_to_file, verify_kmers, verify_kmers_2};
use ahash::RandomState;
use std::path::PathBuf;

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

        let kmers_on_graph : Vec<Kmer> = generate_kmers(graph,kmer_length as u64, Some(max_degree));
        println!("kmers_on_graph: {:#?}", kmers_on_graph);

        // Print to fasta files
        print_seq_to_file(&forward, &PathBuf::from("./output/ref.fa"));
        print_kmers_to_file(&kmers_on_graph, &PathBuf::from("./output/kmers.fa"));
        //verify_kmers(&forward, &kmers_on_graph, &PathBuf::from("./output/alignment.txt"));

        let kmers_on_seq_fwd : Vec<KmerPos> = generate_pos_on_forward(&kmers_on_graph, &forward, &seq_bv, &node_ref);
        //println!("kmers_on_seq_fwd: {:#?}", kmers_on_seq_fwd);
        verify_kmers_2(&forward, &kmers_on_graph, &kmers_on_seq_fwd, &PathBuf::from("./output/kmers_on_fwd_2.txt"));

        // TODO: seed should be fixed?
        let hash_build = RandomState::new();
        let hashes = generate_hash(&kmers_on_graph, &hash_build);
        let kmer_mphf_table = generate_mphf(&hashes);
        println!("{:#?}", kmer_mphf_table);

        index
    }
}