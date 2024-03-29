use std::path::PathBuf;

use clap::ArgMatches;
use gfa::{gfa::GFA, parser::GFAParser};
use handlegraph::hashgraph::HashGraph;

use crate::index::Index;
use log::{info, warn};
use std::env;

pub fn index_main(global_matches: &ArgMatches) {
    let matches = global_matches.subcommand_matches("index").unwrap();

    let in_path_file = matches.value_of("input").unwrap();

    let out_prefix = matches
        .value_of("out-prefix")
        // Keep the same path but remove ".gfa"
        .unwrap_or(&in_path_file[0..in_path_file.len() - 4]);

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

    let sampling_rate = match matches.value_of("sampling-rate") {
        Some(rate) => Some(rate.parse::<u64>().unwrap()),
        _ => None,
    };

    let generate_mappings = matches.is_present("generate-mappings");

    let mappings_path: Option<&str> = matches.value_of("mappings-path");

    let n_threads = matches
        .value_of("n-threads")
        .unwrap_or(&"0") // Use all available threads
        .parse::<usize>()
        .unwrap();

    // Initialize the RUST_LOG env variable for logging.
    // Note: the log.features attribute in Cargo.toml defines what kind of events
    // to show or hide in debug/release mode.
    if env::var("RUST_LOG").is_err() {
        env::set_var("RUST_LOG", "info")
    }
    env_logger::init();
    info!("Start logging");

    /*
    // Create threadpool with the number of threads specified by the user
    match rayon::ThreadPoolBuilder::new().num_threads(n_threads).build_global() {
        Ok(_) => info!("Building index using {} threads!", rayon::current_num_threads()),
        Err(e) => panic!("{}",e)
    };
     */

    // Create HashGraph from GFA
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser.parse_file(&PathBuf::from(in_path_file)).unwrap();
    let graph = HashGraph::from_gfa(&gfa);

    // Build the index for the input graph
    Index::build(
        &graph,
        kmer_length,
        max_furcations,
        max_degree,
        Some(out_prefix),
        sampling_rate,
        generate_mappings,
        mappings_path,
    );
}
