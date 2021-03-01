use clap::{App, Arg};
use gfa::{parser::GFAParser, gfa::GFA};
use std::path::PathBuf;

use handlegraph::hashgraph::HashGraph;

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
        .get_matches();

    let in_path_file = matches
        .value_of("input")
        .expect("Could not parse argument --input");

    let parser = GFAParser::new();
    let gfa : GFA<usize, ()> = parser.parse_file(&PathBuf::from(in_path_file)).unwrap();
    let graph = HashGraph::from_gfa(&gfa);

    println!("{:#?}",graph);

}
