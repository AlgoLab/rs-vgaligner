use clap::ArgMatches;
use gfa::{gfa::GFA, parser::GFAParser};
use handlegraph::hashgraph::HashGraph;
use std::path::PathBuf;
use crate::index::Index;

pub fn index_main(global_matches : &ArgMatches) {
    let matches = global_matches.subcommand_matches("index").unwrap();

    let in_path_file = matches
        .value_of("input")
        .unwrap();

    let out_prefix = matches
        .value_of("out-prefix")
        // Keep the same path but remove ".gfa"
        .unwrap_or(&in_path_file[0..in_path_file.len()-4]);

    let kmer_length = matches
        .value_of("kmer-length")
        .unwrap()
        .parse::<u64>()
        .unwrap();

    let max_furcations = matches
        .value_of("max-furcations")
        .unwrap_or(&"100")
        .parse::<u64>()
        .unwrap();

    let max_degree = matches
        .value_of("max-degree")
        .unwrap_or(&"100")
        .parse::<u64>()
        .unwrap();

    let n_threads = matches
        .value_of("n-threads")
        .unwrap_or(&"0")    // Use all available threads
        .parse::<usize>()
        .unwrap();

    // Create HashGraph from GFA
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser
        .parse_file(&PathBuf::from(in_path_file))
        .unwrap();
    let graph = HashGraph::from_gfa(&gfa);

    // Build the index for the input graph
    Index::build(
        &graph,
        kmer_length,
        max_furcations,
        max_degree,
        Some(out_prefix),
        n_threads,
    );
}