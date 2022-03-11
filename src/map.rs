use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

//use rayon::iter::ParallelIterator;
//use rayon::prelude::{IntoParallelIterator, IntoParallelRefIterator};

use crate::align::{best_alignment_for_query, GAFAlignment};
use crate::chain::{anchors_for_query, chain_anchors, Anchor, Chain};
use crate::index::Index;
use crate::io::QuerySequence;

use log::{info, warn};
use std::time::{Duration, Instant};
use gfa::gfa::GFA;
use gfa::parser::GFAParser;
use handlegraph::hashgraph::HashGraph;
use crate::validate::{create_validation_records, ValidationRecord, write_validation_to_file};

/// Map the [input] reads against the [index], in order to obtain Chains
/// (and eventually POA Alignments is [also_align]), which can be printed to console
/// with [out-console] or store the in a file with prefix [out_prefix].
/// Additional parameters include:
/// - [bandwidth] which specifies the search distance
/// - [max_gap] which specifies the maximum distance of anchors that can be chained
/// - [chain_min_n_anchors] which specifies the minimum number of anchors that make up a chain
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
    align_best_n: u64,
) {
    info!("Found {} reads!", inputs.len());

    let start_chaining = Instant::now();
    // Collect chains obtained from each input sequence
    let chains: Vec<Vec<Chain>> = inputs
        .into_iter()
        .map(|query| {
            // First find the anchors, aka exact matches between
            // seqs and kmers in the index
            let mut seq_anchors: Vec<Anchor> = anchors_for_query(index, query);

            /*
            println!("noderef: {:#?}", index.node_ref);
            for anchor in &seq_anchors {
                println!("Anchor: {}, start_id: {}, end_id: {}, start_pos: {:#?}, end_pos: {:#?}",
                         anchor.id,
                         index.handle_from_seqpos(&anchor.target_begin).id(),
                         index.handle_from_seqpos(&anchor.target_end).id(),
                         anchor.target_begin,
                         anchor.target_end
                );
            }
             */

            // Chain close anchors together to find longer matches
            let mut seq_chains: Vec<Chain> = chain_anchors(
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

            //println!("Chains: {:#?}", seq_chains);
            /*
            for chain in &seq_chains {
                for anchor in &chain.anchors {
                    println!("Anchor: {}, \
                         start_id: {}, end_id: {}, \
                         start_pos: {:#?}, end_pos: {:#?}",
                             anchor.id,
                             index.handle_from_seqpos(&anchor.target_begin).id(),
                             index.handle_from_seqpos(&anchor.target_end).id(),
                             anchor.target_begin,
                             anchor.target_end
                    );
                }
                println!("\n\n")
            }
             */

            // Sort the chains according to their number of anchors
            // TODO: using the score would be better...
            seq_chains.sort_by(|a, b| a.anchors.len().cmp(&b.anchors.len()));

            seq_chains
        })
        .collect();
    info!("Chaining took: {} ms", start_chaining.elapsed().as_millis());

    // Store chains as GAF Strings, and store them to a output GAF file,
    // and optionally to console
    info!(
        "Found {} chains!",
        chains.iter().map(|x| x.len()).sum::<usize>()
    );

    //info!("Chains are: {:#?}", chains);

    let chains_gaf: Vec<GAFAlignment> = chains
        .iter()
        .flat_map(|query_chains| {
            query_chains.iter().map(|c|
                    // TODO: overflow when creating alignment from chain
                    match c.is_placeholder {
                        false => GAFAlignment::from_chain(c, index),
                        true => GAFAlignment::from_placeholder_chain(c)
                })
        })
        .collect();

    if let Some(prefix) = out_prefix {
        let file_name = match prefix.ends_with(".gaf") {
            true => prefix.to_string(),
            false => prefix.clone().to_string() + "-chains" + ".gaf",
        };

        match write_gaf_to_file(&chains_gaf, file_name.to_string()) {
            Err(e) => panic!("{}", e),
            _ => info!("Chains stored correctly in {}!", file_name),
        }
    }

    if write_console {
        for gaf_str in chains_gaf {
            info!("{:#?}", gaf_str);
        }
    }

    // Perform the alignment with rs-abPOA, using the chains as a reference
    if also_align {
        let start_alignment = Instant::now();
        let alignments: Vec<GAFAlignment> = chains
            .iter()
            .map(|query_chains| best_alignment_for_query(index, query_chains, align_best_n))
            .collect();
        info!(
            "Alignment took: {} ms",
            start_alignment.elapsed().as_millis()
        );

        info!("Found {} alignments!", alignments.len());
        if let Some(prefix) = out_prefix {
            let file_name = match prefix.ends_with(".gaf") {
                true => prefix.to_string(),
                false => prefix.clone().to_string() + "-alignments" + ".gaf",
            };

            match write_gaf_to_file(&alignments, file_name.to_string()) {
                Err(e) => panic!("{}", e),
                _ => info!("Alignments stored correctly in {}!", file_name),
            }
        }

        /*
        // TODO: add as parameter
        let also_validate = true;
        let input_file = "";
        let validation_file_name = "";

        if also_validate {
            // Create HashGraph from GFA
            let parser = GFAParser::new();
            let gfa: GFA<usize, ()> = parser
                .parse_file(&PathBuf::from(input_file))
                .unwrap();
            let graph = HashGraph::from_gfa(&gfa);

            // Obtain validation records
            let val_records : Vec<ValidationRecord> = create_validation_records(&graph,&alignments, &inputs);

            // Write validation records to file
            match write_validation_to_file(&val_records, validation_file_name.to_string()) {
                Err(e) => panic!("{}", e),
                _ => info!("Alignments stored correctly in {}!", validation_file_name),
            }
        }
         */

        if write_console {
            for gaf_str in alignments {
                println!("{:#?}", gaf_str);
            }
        }
    }
}

/// Store [chains_gaf_strings] in the file with name [file_name]
fn write_gaf_to_file(gaf_alignments: &Vec<GAFAlignment>, file_name: String) -> std::io::Result<()> {
    let gaf_strings: Vec<String> = gaf_alignments
        .iter()
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
    use std::path::PathBuf;

    use gfa::gfa::GFA;
    use gfa::parser::GFAParser;
    use handlegraph::hashgraph::HashGraph;

    use crate::index::Index;
    use crate::io::read_seqs_from_file;
    use crate::map::map_reads;

    #[test]
    fn test_map_no_alignment() {
        let parser = GFAParser::new();
        let gfa: GFA<usize, ()> = parser
            .parse_file(&PathBuf::from("./test/test.gfa"))
            .unwrap();
        let graph = HashGraph::from_gfa(&gfa);

        let index = Index::build(&graph, 11, 100, 100, None, None, false, None);
        let query = read_seqs_from_file("./test/single-read-test.fa").unwrap();

        map_reads(
            &index, &query, 50, 1000, 3, 0.5f64, 0.1, 60.0f64, false, None, false, 100,
        );
    }

    #[test]
    fn test_map_with_alignment() {
        let parser = GFAParser::new();
        let gfa: GFA<usize, ()> = parser
            .parse_file(&PathBuf::from("./test/test.gfa"))
            .unwrap();
        let graph = HashGraph::from_gfa(&gfa);

        let index = Index::build(&graph, 11, 100, 100, None, None, false, None);
        let query = read_seqs_from_file("./test/single-read-test.fa").unwrap();

        map_reads(
            &index, &query, 50, 1000, 3, 0.5f64, 0.1, 60.0f64, false, None, true, 100,
        );
    }
}
