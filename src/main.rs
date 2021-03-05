use std::path::PathBuf;
use gfa::{parser::GFAParser, gfa::GFA};
use handlegraph::handlegraph::HandleGraph;
use handlegraph::hashgraph::HashGraph;
use clap::{App, Arg};

mod index;
mod utils;
mod dna;

//use crate::index::{Index,IndexUtilities};
use crate::index::Index;
use crate::dna::reverse_complement;

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

    // STEP 1: Build the index for the input graph
    let graph_index = Index::build(&graph, &kmer_length, &100, &100, &1.0, &"./output/file");

}
