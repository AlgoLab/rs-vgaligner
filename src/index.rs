use bv::BitVec;
use handlegraph::handle:: {Handle,Edge,Direction};
use handlegraph::hashgraph::HashGraph;
use handlegraph::handlegraph::HandleGraph;
use bstr::ByteSlice;

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

pub struct node_ref {
    seq_idx : u64,
    edge_idx : u64,
    count_prev : u64,
}

// TODO: add remaining traits etc.

pub trait IndexUtilities {
    fn build(graph : &HashGraph,
             kmer_length : &u64,
             max_furcations : &u64,
             max_degree : &u64,
             sampling_rate : &f32,
             out_prefix : &str) -> Index;
    //fn load() -> Index;
}

impl IndexUtilities for Index {
    fn build(graph : &HashGraph, kmer_length : &u64, max_furcations : &u64, max_degree : &u64, sampling_rate : &f32, out_prefix : &str) -> Index {
        let number_nodes = graph.graph.len();

        // Find total length of the sequences in the graph
        let mut total_length = 0;
        for value in graph.handles_iter() {
            let node = graph.get_node(&value.id()).unwrap();
            total_length += node.sequence.len() as u64;
        }

        let mut seq_bw : BitVec = BitVec::new_fill(false,total_length+1);

        let mut seq_idx : u64 = 0;

        let mut n_edges = 0; //TODO: review

        graph.handles_iter().for_each(|h| {
            seq_bw.set(seq_idx, true);

            let node = graph.get_node(&h.id()).unwrap();
            let seq = node.sequence.clone().to_string();

            let mut reference = node_ref {
                seq_idx,
                edge_idx : n_edges,
                count_prev : 0
            };

            graph.handle_edges_iter(h, Direction::Right).for_each(|p|{
                reference.count_prev += 1;
            });

            n_edges += reference.count_prev;

            graph.handle_edges_iter(h, Direction::Left).for_each(|p|{
                n_edges += 1;
            });

            seq_idx += seq.len() as u64;
        });

        Index {
            kmer_length: *kmer_length,
            sampling_length: 0,
            seq_length: total_length,
            seq_fwd: Vec::new(),
            seq_rev: Vec::new(),
            seq_bv: seq_bw,
            seq_by_rank: BitVec::new(),
            n_edges: n_edges,
            edges: Vec::new(),
            n_nodes: number_nodes as u64,
            nodes: Vec::new(),
            n_kmers: 0,
            n_kmer_pos: 0,
            kmer_pos_ref: Vec::new(),
            kmer_pos_table: Vec::new(),
            loaded: true
        }
    }
}