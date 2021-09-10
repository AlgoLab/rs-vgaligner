use clap::ArgMatches;
use crate::map::map_reads;
use crate::index::Index;
use crate::io::read_seqs_from_file;

pub fn map_main(global_matches : &ArgMatches) {
    let matches = global_matches.subcommand_matches("map").unwrap();

    let idx_prefix = matches
        .value_of("input")
        .unwrap();

    let in_path_file = matches
        .value_of("input-file")
        .unwrap();

    let out_prefix = matches
        .value_of("out-prefix")
        // Keep the same path but remove ".gfa"
        .unwrap_or(&idx_prefix);

    let max_gap_length = matches
        .value_of("max-gap-length")
        .unwrap_or(&"1000")
        .parse::<u64>()
        .unwrap();

    let max_mismatch_rate = matches
        .value_of("max-mismatch-rate")
        .unwrap_or(&"0.1")
        .parse::<f64>()
        .unwrap();

    let chain_min_n_anchors = matches
        .value_of("chain-min-anchors")
        .unwrap_or(&"3")
        .parse::<u64>()
        .unwrap();

    // TODO: review
    let _align_best_n = matches
        .value_of("align-best-n")
        .unwrap_or(&"100")
        .parse::<u64>()
        .unwrap();

    let write_chains = matches
        .is_present("write-chains");

    let dont_align = matches
        .is_present("dont-align");

    let n_threads = matches
        .value_of("n-threads")
        .unwrap_or(&"0")    // Use all available threads
        .parse::<usize>()
        .unwrap();

    let index = Index::load_from_prefix(idx_prefix.to_string());

    let query = read_seqs_from_file(&in_path_file).unwrap();

    map_reads(&index, &query, 50, max_gap_length,
              chain_min_n_anchors, 0.5f64,
              max_mismatch_rate, 60.0f64,
              write_chains, Some(out_prefix), dont_align,
              n_threads
    );
}