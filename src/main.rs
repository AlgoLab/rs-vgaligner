use clap::{App, Arg};
use gfa::{gfa::GFA, parser::GFAParser};
use handlegraph::hashgraph::HashGraph;
use std::path::PathBuf;

pub mod dna;
pub mod index;
pub mod io;
pub mod kmer;
pub mod utils;
mod serialization;

use crate::index::Index;

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
            Arg::with_name("out-prefix")
                .short("o")
                .long("out")
                .value_name("STRING")
                .help("Save our index with this prefix (defaults to input file name and path)")
                .required(false)
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
        .arg(
            Arg::with_name("max-furcations")
                .short("e")
                .long("max-furcations")
                .value_name("INTEGER")
                .help("Break at edges that would induce this many furcations in a kmer")
                .required(false)
                .takes_value(true),
        )
        .arg(
            Arg::with_name("max-degree")
                .short("D")
                .long("max-degree")
                .value_name("INTEGER")
                .help("Remove nodes that have degree greater that this level")
                .required(false)
                .takes_value(true),
        )
        .arg(
            Arg::with_name("sampling-rate")
                .short("r")
                .long("sampling-rate")
                .value_name("FLOAT")
                .help("Build a modimizer index by keeping kmers whose hash % round(RATE/1) == 0")
                .required(false)
                .takes_value(true),
        )
        .get_matches();

    let in_path_file = matches.value_of("input").unwrap();

    let out_prefix = matches.value_of("out-prefix").unwrap_or_else(|| &"");

    let kmer_length = matches
        .value_of("kmer-length")
        .unwrap()
        .parse::<u64>()
        .unwrap();

    let max_furcations = matches
        .value_of("max-furcations")
        .unwrap_or_else(|| &"100")
        .parse::<u64>()
        .unwrap();

    let max_degree = matches
        .value_of("max-degree")
        .unwrap_or_else(|| &"100")
        .parse::<u64>()
        .unwrap();

    let sampling_rate = matches
        .value_of("sampling-rate")
        .unwrap_or_else(|| &"1.0")
        .parse::<f32>()
        .unwrap();

    // Create HashGraph from GFA
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser.parse_file(&PathBuf::from(in_path_file)).unwrap();
    let graph = HashGraph::from_gfa(&gfa);

    // STEP 1: Build the index for the input graph
    let graph_index = Index::build(
        &graph,
        kmer_length,
        max_furcations,
        max_degree,
        sampling_rate,
        out_prefix,
    );
}
