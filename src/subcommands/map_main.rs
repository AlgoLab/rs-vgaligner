use clap::ArgMatches;
use crate::map::map_reads;
use crate::index::Index;
use crate::io::read_seqs_from_file;

pub fn map_main(global_matches : &ArgMatches) {
    let matches = global_matches.subcommand_matches("map").unwrap();

    let idx_prefix = matches.value_of("input").unwrap();

    let in_path_file = matches.value_of("input-file").unwrap();

    let out_prefix = matches.value_of("out-prefix");

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

    let _chain_overlap_max = matches
        .value_of("chain-overlap-max")
        .unwrap_or(&"0.1")
        .parse::<f32>()
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

    let _write_superchains = matches
        .is_present("write-superchains");

    let dont_align = matches
        .is_present("dont-align");

    let index = Index::load_from_prefix(idx_prefix.to_string());

    let query = read_seqs_from_file(&in_path_file).unwrap();

    map_reads(&index, &query,50, max_gap_length,
              chain_min_n_anchors, 0.5f64,
              max_mismatch_rate, 60.0f64,
              write_chains, out_prefix, dont_align);
}