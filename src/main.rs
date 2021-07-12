#[macro_use]
extern crate clap;
use clap::App;

mod subcommands {
    pub mod index_main;
    pub mod map_main;
}

pub mod dna;
pub mod index;
pub mod map;
pub mod io;
pub mod kmer;
mod serialization;
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