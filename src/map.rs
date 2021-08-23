use crate::index::Index;
use crate::io::InputSequence;
use crate::chain::{Anchor, Chain, write_chain_gaf, anchors_for_query, chain_anchors, extract_graph_po_range};
use std::fs::File;
use std::io::Write;
use rayon::iter::{IntoParallelRefIterator, IntoParallelIterator};

/// Map the [input] reads against the [index].
// TODO: add explaination to other parameters
pub fn map_reads(index : &Index, inputs : &Vec<InputSequence>,
                 bandwidth : u64, max_gap : u64,
                 chain_min_n_anchors : u64, secondary_chain_threshold : f64,
                 max_mismatch_rate : f64, max_mapq : f64,
                 write_chains : bool, out_prefix: Option<&str>) {

    let mut chains : Vec<Chain> = Vec::new();
    let mut chains_gaf_strings : Vec<String> = Vec::new();

    // TODO: use start_id to disambiguate anchors in different for iterations?
    // i.e. the anchor 23 appeared multiple times in chains, that was because
    // in each for iteration I use as anchor_ids (0..seq_anchors.len()), therefore
    // the id is only relative to that specific seq_anchors
    //let mut start_id : anchor_id = 0;

    for query in inputs {
        // First find the anchors, aka exact matches between
        // seqs and kmers in the index
        let mut seq_anchors: Vec<Anchor> = anchors_for_query(index, query);
        
        // Chain close anchors together to find longer matches
        let seq_chains : Vec<Chain> = chain_anchors(&mut seq_anchors, index.kmer_length,bandwidth,
                                         max_gap, chain_min_n_anchors, secondary_chain_threshold,
                                         max_mismatch_rate, max_mapq);

        // Convert chains to strings and do the same
        let curr_chains_gaf_strings : Vec<String> = seq_chains
            .iter()
            .map(|chain| write_chain_gaf(chain, index, &query.name, query.seq.len()))
            .collect();

        // Find partially ordered ranges implied by chains
        extract_graph_po_range(index, &seq_chains, query);

        // Add results to the ones from previous iterations
        chains.extend(seq_chains.into_iter());
        chains_gaf_strings.extend(curr_chains_gaf_strings.into_iter());

    }

    if write_chains {
        if chains.is_empty() {
            println!("No chains found!");
        } else {
            match out_prefix {
                Some(prefix) => {
                    match write_chains_to_file(&chains_gaf_strings, prefix.clone().to_string() + ".gaf") {
                        Err(e) => panic!("{}",e),
                        _ => println!("Chains stored correctly!"),
                    }
                },
                _ => {
                    for gaf_str in chains_gaf_strings {
                        println!("{}", gaf_str);
                    }
                }
            }
        }

    }

}

/// Store [chains_gaf_strings] in the file with name [file_name]
fn write_chains_to_file(
    chains_gaf_strings : &Vec<String>,
    file_name: String,
) -> std::io::Result<()> {
    let mut file = File::create(&file_name).expect(&format!("Couldn't create file {}",&file_name));
    file.write_all(&chains_gaf_strings.join("").as_bytes()).expect(&format!("Couldn't write to file {}", &file_name));
    Ok(())
}

mod test {
    use super::*;
    use crate::io::read_seqs_from_file;
    use gfa::parser::GFAParser;
    use gfa::gfa::GFA;
    use handlegraph::hashgraph::HashGraph;
    use std::path::PathBuf;

    #[test]
    fn test_map() {
        let parser = GFAParser::new();
        let gfa: GFA<usize, ()> = parser.parse_file(&PathBuf::from("./test/test.gfa")).unwrap();
        let graph = HashGraph::from_gfa(&gfa);

        let index = Index::build(&graph, 11, 100,
                                 100, 7.0, None);
        let query = read_seqs_from_file("./test/test.fa").unwrap();

        map_reads(&index, &query,50, 1000,
                  3, 0.5f64,
                  0.1, 60.0f64,
                  false, None);
    }
}