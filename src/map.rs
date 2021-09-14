use crate::chain::{anchors_for_query, chain_anchors, write_chain_gaf, Anchor, Chain};
use crate::index::Index;
use crate::io::QuerySequence;

use handlegraph::hashgraph::HashGraph;
use rayon::prelude::{IntoParallelIterator, IntoParallelRefIterator};
use std::fs::File;
use std::io::Write;

use crate::align::{best_alignment_for_query, obtain_base_level_alignment, GAFAlignment};
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
    write_console: bool,
    out_prefix: Option<&str>,
    also_align: bool,
) {
    println!("Found {} reads!", inputs.len());

    // Collect chains obtained from each input sequence
    let chains: Vec<Vec<Chain>> = inputs
        .into_par_iter()
        .filter_map(|query| {
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

            if !seq_chains.is_empty() {
                Some(seq_chains)
            } else {
                None
            }
        })
        .collect();

    // Store chains as GAF Strings, and store them to a output GAF file,
    // and optionally to console
    if chains.is_empty() {
        println!("No chain found!");
    } else {
        println!(
            "Found {} chains!",
            chains.par_iter().map(|x| x.len()).sum::<usize>()
        );

        let chains_gaf: Vec<GAFAlignment> = chains
            .par_iter()
            .flat_map(|query_chains| query_chains.par_iter().map(|c| GAFAlignment::from_chain(c)))
            .collect();

        if let Some(prefix) = out_prefix {
            match write_gaf_to_file(&chains_gaf, prefix.clone().to_string() + "-chains" + ".gaf") {
                Err(e) => panic!("{}", e),
                _ => println!(
                    "Chains stored correctly in {}-chains.gaf!",
                    prefix.clone().to_string()
                ),
            }
        }

        if write_console {
            for gaf_str in chains_gaf {
                println!("{:#?}", gaf_str);
            }
        }
    }

    // Perform the alignment with rs-abPOA, using the chains as a reference
    if also_align {
        let alignments: Vec<GAFAlignment> = chains
            .par_iter()
            .map(|query_chains| best_alignment_for_query(index, query_chains))
            .collect();

        if alignments.is_empty() {
            println!("No alignment found!");
        } else {
            println!("Found {} alignments!", alignments.len());
            if let Some(prefix) = out_prefix {
                match write_gaf_to_file(
                    &alignments,
                    prefix.clone().to_string() + "-alignments" + ".gaf",
                ) {
                    Err(e) => panic!("{}", e),
                    _ => println!(
                        "Alignments stored correctly in {}-alignments.gaf!",
                        prefix.clone().to_string()
                    ),
                }
            }

            if write_console {
                for gaf_str in alignments {
                    println!("{:#?}", gaf_str);
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
    fn test_map_no_alignment() {
        let parser = GFAParser::new();
        let gfa: GFA<usize, ()> = parser
            .parse_file(&PathBuf::from("./test/test.gfa"))
            .unwrap();
        let graph = HashGraph::from_gfa(&gfa);

        let index = Index::build(&graph, 11, 100, 100, None);
        let query = read_seqs_from_file("./test/single-read-test.fa").unwrap();

        map_reads(
            &index, &query, 50, 1000, 3, 0.5f64, 0.1, 60.0f64, false, None, false,
        );
    }

    #[test]
    fn test_map_with_alignment() {
        let parser = GFAParser::new();
        let gfa: GFA<usize, ()> = parser
            .parse_file(&PathBuf::from("./test/test.gfa"))
            .unwrap();
        let graph = HashGraph::from_gfa(&gfa);

        let index = Index::build(&graph, 11, 100, 100, None);
        let query = read_seqs_from_file("./test/single-read-test.fa").unwrap();

        map_reads(
            &index, &query, 50, 1000, 3, 0.5f64, 0.1, 60.0f64, false, None, true,
        );
    }
}
