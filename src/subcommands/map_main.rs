use clap::ArgMatches;

use crate::index::Index;
use crate::io::read_seqs_from_file;
use crate::map::map_reads;

pub fn map_main(global_matches : &ArgMatches) {
    let matches = global_matches.subcommand_matches("map").unwrap();

    let idx_prefix = matches
        .value_of("index")
        .unwrap();

    let in_path_file = matches
        .value_of("input-file")
        .unwrap();

    let out_prefix = matches
        .value_of("out-prefix")
        // Keep the same path but remove ".fa/fasta/fq/fastq"
        .unwrap_or_else(||
            if in_path_file.ends_with("fa") || in_path_file.ends_with("fasta") {
                &in_path_file[0..in_path_file.len()-3]
            } else {
                &in_path_file[0..in_path_file.len()-4]
            });

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

    let align_best_n = matches
        .value_of("align-best-n")
        .unwrap_or(&"1")
        .parse::<u64>()
        .unwrap();

    let write_console = matches
        .is_present("write-console");

    let also_align = matches
        .is_present("also-align");

    let n_threads = matches
        .value_of("n-threads")
        .unwrap_or(&"0")    // Use all available threads
        .parse::<usize>()
        .unwrap();

    // Create threadpool with the number of threads specified by the user
    match rayon::ThreadPoolBuilder::new().num_threads(n_threads).build_global() {
        Ok(_) => println!("Mapping reads to Index using {} threads!", rayon::current_num_threads()),
        Err(e) => panic!("{}",e)
    };

    let index: Index = match idx_prefix.ends_with(".idx") {
        true => Index::load_from_file(idx_prefix.to_string()),
        false => Index::load_from_prefix(idx_prefix.to_string()),
    };

    let query = read_seqs_from_file(&in_path_file).unwrap();

    map_reads(&index, &query, 50, max_gap_length,
              chain_min_n_anchors, 0.5f64,
              max_mismatch_rate, 60.0f64,
              write_console, Some(out_prefix), also_align, align_best_n);
}