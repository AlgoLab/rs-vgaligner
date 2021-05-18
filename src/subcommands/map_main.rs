use clap::ArgMatches;
use crate::map::map_reads;

pub fn map_main(global_matches : &ArgMatches) {
    let matches = global_matches.subcommand_matches("map").unwrap();
    map_reads();
}