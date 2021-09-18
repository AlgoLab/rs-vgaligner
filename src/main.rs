#[macro_use]
extern crate clap;

use clap::App;

// Subcommands
mod subcommands {
    pub mod index_main;
    pub mod map_main;
}

// Main program logic
pub mod align;
pub mod index;
pub mod map;

// Main data structures used
pub mod chain;
pub mod kmer;

// Io
pub mod io;
mod serialization;

// Utils
pub mod dna;
pub mod utils;

fn main() {
    let yaml = load_yaml!("./subcommands/cli.yml");
    let matches = App::from_yaml(yaml).get_matches();

    match matches.subcommand_name() {
        Some("index") => subcommands::index_main::index_main(&matches),
        Some("map") => subcommands::map_main::map_main(&matches),
        _ => println!("Missing subcommand, please add [index|map]"),
    }
}
