use bv::BitVec;
use handlegraph::handle::{Handle,Edge,Direction,NodeId};
use handlegraph::hashgraph::HashGraph;
use handlegraph::handlegraph::HandleGraph;
use bstr::ByteSlice;
use std::fs::File;
use std::io::{Write, Read};
use crate::utils::{get_bv_rank, find_sequence};
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

        // TODO: add node_ref
        let forward = find_sequence(graph, &mut seq_bv);
        let reverse = reverse_complement(&forward.as_str());

        println!("Forward is: {}", forward);
        println!("Reverse is: {}", reverse);
        println!("BV is: {:#?}", seq_bv);

        let kmers_on_graph = generate_kmers(graph,kmer_length as u64, Some(max_degree));
        //let kmers_on_seq_fwd : Vec<KmerSeq> = generate_pos_on_fwd(kmers_on_graph, seq_bv, node_ref);
        //let hashes = generate_kmers_hash(&kmers_on_seq_fwd);

        println!("kmers_on_graph: {:#?}", kmers_on_graph);
        //println!("hashes: {:#?}", hashes);

        //for val in hashes.keys() {
        //    println!("{}",val);
        //}


        /*
        // Create files

        // This file will contain the forward sequences of the graph
        let seq_fwd_filename : String = out_prefix.to_owned() + ".sqf";
        let mut seq_fwd_f = File::create(seq_fwd_filename).unwrap();

        // This file will contain the reverse sequences of the graph
        let seq_rev_filename : String = out_prefix.to_owned() + ".sqr";
        let mut seq_rev_f = File::create(seq_rev_filename).unwrap();

        // This file will contain the incoming edges of a given node
        let edge_filename : String = out_prefix.to_owned() + ".gye";
        let mut edge_f = File::create(edge_filename).unwrap();

        //This file will contain the outgoing edges of a given node
        let node_ref_filename : String = out_prefix.to_owned() + ".gyn";
        let mut node_ref_f = File::create(node_ref_filename).unwrap();


        // Iterate over the graph to fill in the files defined above
        let mut seq_idx : u64 = 0;
        let mut n_edges = 0; //TODO: review -- not initialized in C++?

        graph.handles_iter().for_each(|h| {

            // Mark the node in the bitvector
            seq_bv.set(seq_idx, true);

            // Get the sequence of the node
            let node = graph.get_node(&h.id()).unwrap();
            let seq = node.sequence.clone().to_string();

            // Write the forward and reverse complemented version
            seq_fwd_f.write(&seq.as_bytes());

            //seq_rev_f.write(&reverse_complement(&seq).as_bytes());
            // Will be reverse complemented later on
            seq_rev_f.write(&seq.as_bytes());

            let mut reference = Node_ref {
                seq_idx,
                edge_idx : n_edges,
                count_prev : 0
            };

            graph.handle_edges_iter(h, Direction::Left).for_each(|p|{
                edge_f.write(p.id().to_string().as_bytes());
                reference.count_prev += 1;
            });

            n_edges += reference.count_prev;
            let reference_string= format!("{},{},{}\n",
                                          reference.seq_idx,
                                          reference.edge_idx,
                                          reference.count_prev);
            node_ref_f.write(reference_string.as_bytes());

            graph.handle_edges_iter(h, Direction::Right).for_each(|p|{
                edge_f.write(p.id().to_string().as_bytes());
                n_edges += 1;
            });

            seq_idx += seq.len() as u64;
        });

        // Write a marker reference, to simplify counting of edges

        let mut reference = Node_ref {
            seq_idx,
            edge_idx : n_edges,
            count_prev : 0
        };

        let reference_string= format!("{},{},{}\n",
                                      reference.seq_idx,
                                      reference.edge_idx,
                                      reference.count_prev);
        node_ref_f.write(reference_string.as_bytes());

        assert_eq!(reference.seq_idx, total_length);


        // Save bw and rank as files
        seq_bv.set(seq_idx, true); //end marker

        //This file will the bitvector
        let seq_bv_filename : String = out_prefix.to_owned() + ".sbv";
        let mut seq_bv_f = File::create(seq_bv_filename).unwrap();

        // Compute bitvector rank
        let mut seq_bv_rank = get_bv_rank(&seq_bv);

        let mut seq_fwd = String::new();

        seq_fwd_f.read_to_string(&mut seq_fwd);
        let seq_rev = reverse_complement(seq_fwd.as_str());

        println!("Forward: {:#?}", seq_fwd);
        println!("Reverse: {:#?}", seq_rev);

        //println!("Seq bv: {:#?}",seq_bv);
        //println!("Seq rank: {:#?}", seq_bv_rank);

        //seq_bv.serialize(&seq_bv_f);
        //seq_by_rank.serialize(&seq_bv_f);
         */

        index
    }
}