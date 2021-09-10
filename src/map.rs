use crate::chain::{anchors_for_query, chain_anchors, write_chain_gaf, Anchor, Chain};
use crate::index::Index;
use crate::io::QuerySequence;

use handlegraph::hashgraph::HashGraph;
use rayon::prelude::{IntoParallelIterator, IntoParallelRefIterator};
use std::fs::File;
use std::io::Write;

use crate::align::{obtain_base_level_alignment, GAFAlignment};
use rayon::iter::ParallelIterator;

/// Map the [input] reads against the [index].
// TODO: add explaination to other parameters
pub fn map_reads(
    index: &Index,
    inputs: &Vec<QuerySequence>,
    bandwidth: u64,
    max_gap: u64,
    chain_min_n_anchors: u64,
    secondary_chain_threshold: f64,
    max_mismatch_rate: f64,
    max_mapq: f64,
    write_chains: bool,
    out_prefix: Option<&str>,
    dont_align: bool,
    n_threads: usize,
) {
    // Create threadpool with the number of threads specified by the user
    match rayon::ThreadPoolBuilder::new().num_threads(n_threads).build_global() {
        Ok(_) => println!("Mapping reads to Index using {} threads", rayon::current_num_threads()),
        Err(_) => println!("Threadpool already initialized!")
    };

    // TODO: use start_id to disambiguate anchors in different for iterations?
    // i.e. the anchor 23 appeared multiple times in chains, that was because
    // in each for iteration I use as anchor_ids (0..seq_anchors.len()), therefore
    // the id is only relative to that specific seq_anchors
    //let mut start_id : anchor_id = 0;

    let chains: Vec<Chain> = inputs
        .into_par_iter()
        .flat_map(|query| {
            // First find the anchors, aka exact matches between
            // seqs and kmers in the index
            let mut seq_anchors: Vec<Anchor> = anchors_for_query(index, query);

            // Chain close anchors together to find longer matches
            let seq_chains: Vec<Chain> = chain_anchors(
                &mut seq_anchors,
                index.kmer_length,
                bandwidth,
                max_gap,
                chain_min_n_anchors,
                secondary_chain_threshold,
                max_mismatch_rate,
                max_mapq,
                query,
            );

            seq_chains
        })
        .collect();

    if write_chains {
        if chains.is_empty() {
            println!("No chain found!");
        } else {
            let chains_gaf: Vec<GAFAlignment> = chains
                .par_iter()
                .map(|c| GAFAlignment::from_chain(c))
                .collect();
            match out_prefix {
                Some(prefix) => {
                    match write_gaf_to_file(&chains_gaf, prefix.clone().to_string() + ".gaf") {
                        Err(e) => panic!("{}", e),
                        _ => println!("Chains stored correctly!"),
                    }
                }
                _ => {
                    for gaf_str in chains_gaf {
                        println!("{:#?}", gaf_str);
                    }
                }
            }
        }
    }

    if !dont_align {
        let alignments: Vec<GAFAlignment> = chains
            .par_iter()
            .map(|chain| obtain_base_level_alignment(index, chain))
            .collect();

        if alignments.is_empty() {
            println!("No alignment found!");
        } else {
            match out_prefix {
                Some(prefix) => {
                    match write_gaf_to_file(
                        &alignments,
                        prefix.clone().to_string() + "-aligned" + ".gaf",
                    ) {
                        Err(e) => panic!("{}", e),
                        _ => println!("Alignments stored correctly!"),
                    }
                }
                _ => {
                    for gaf_str in alignments {
                        println!("{:#?}", gaf_str);
                    }
                }
            }
        }
    }
}

/// Store [chains_gaf_strings] in the file with name [file_name]
fn write_gaf_to_file(gaf_alignments: &Vec<GAFAlignment>, file_name: String) -> std::io::Result<()> {
    let gaf_strings: Vec<String> = gaf_alignments
        .par_iter()
        .map(|aln| aln.to_string())
        .collect();
    let mut file =
        File::create(&file_name).unwrap_or_else(|_| panic!("Couldn't create file {}", &file_name));
    file.write_all(&gaf_strings.join("").as_bytes())
        .unwrap_or_else(|_| panic!("Couldn't write to file {}", &file_name));
    Ok(())
}

#[cfg(test)]
mod test {
    use crate::index::Index;
    use crate::io::read_seqs_from_file;
    use crate::map::map_reads;
    use gfa::gfa::GFA;
    use gfa::parser::GFAParser;
    use handlegraph::hashgraph::HashGraph;
    use std::path::PathBuf;

    #[test]
    fn test_map() {
        let parser = GFAParser::new();
        let gfa: GFA<usize, ()> = parser
            .parse_file(&PathBuf::from("./test/test.gfa"))
            .unwrap();
        let graph = HashGraph::from_gfa(&gfa);

        let index = Index::build(&graph, 11, 100, 100, None, 0);
        let query = read_seqs_from_file("./test/single-read-test.fa").unwrap();

        map_reads(
            &index, &query, 50, 1000, 3, 0.5f64, 0.1, 60.0f64, false, None, true, 0,
        );
    }
}
