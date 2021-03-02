use std::path::PathBuf;
use gfa::{parser::GFAParser, gfa::GFA};
use handlegraph::handlegraph::HandleGraph;
use handlegraph::hashgraph::HashGraph;
use clap::{App, Arg};
use bv::BitVec;
use bstr::ByteSlice;

pub struct node_ref {
    seq_idx : u64, // index among sequences
    edge_idx : u64, // index among edges
    count_prev : u64,
}

fn main() {

    let matches = App::new("rs-vgalign")
        .version("0.1")
        .author("Francesco Porto <francesco.porto97@gmail.com>")
        .about("Aligns reads to a Variation Graph")
        .arg(
            Arg::with_name("input")
                .short("i")
                .long("input")
                .value_name("FILE")
                .help("Sets the input file to use")
                .required(true)
                .takes_value(true),
        )
        .arg(
            Arg::with_name("kmer-length")
                .short("k")
                .long("kmer-length")
                .value_name("INTEGER")
                .help("Sets the size of each k-mer")
                .required(true)
                .takes_value(true),
        )
        .get_matches();

    let in_path_file = matches
        .value_of("input")
        .unwrap();

    let kmer_length = matches
        .value_of("kmer-length")
        .unwrap()
        .parse::<u64>()
        .unwrap();

    // Create HashGraph from GFA
    let parser = GFAParser::new();
    let gfa : GFA<usize, ()> = parser.parse_file(&PathBuf::from(in_path_file)).unwrap();
    let graph = HashGraph::from_gfa(&gfa);

    let number_nodes = graph.graph.len();

    // Find total length of the sequences in the graph
    let mut total_length = 0;
    for value in graph.handles_iter() {
        let node = graph.get_node(&value.id()).unwrap();
        total_length += node.sequence.len() as u64;
    }

    let mut seq_bw : BitVec = BitVec::new_fill(false,total_length+1);

    let mut seq_idx : u64 = 0;

    graph.handles_iter().map(|h| {
        seq_bw[seq_idx] = true;

        let node = graph.get_node(&h.id()).unwrap();
        let seq = node.sequence.clone().to_string();

        let reference = node_ref {
            seq_idx : seq_idx,
            edge_idx : n_edges,
            count_prev : 0
        };

        graph.edges_iter().map(|e| {

        });

        n_edges += reference.count_prev;

        graph.edges_iter().map(|e| {

        });

        seq_idx += seq.len();
    });


    println!("{:#?}",bitvector);

}
