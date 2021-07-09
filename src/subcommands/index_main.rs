use clap::ArgMatches;
use gfa::{gfa::GFA, parser::GFAParser};
use handlegraph::hashgraph::HashGraph;
use std::path::PathBuf;
use crate::index::Index;
use handlegraph::handlegraph::HandleGraph;

pub fn index_main(global_matches : &ArgMatches) {
    let matches = global_matches.subcommand_matches("index").unwrap();

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
    let _graph_index = Index::build(
        &graph,
        kmer_length,
        max_furcations,
        max_degree,
        sampling_rate,
        out_prefix.to_string(),
    );
}