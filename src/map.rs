use crate::index::Index;
use crate::io::InputSequence;
use crate::chain::{Anchor, Chain, write_chain_gaf, anchors_for_query, chain_anchors};
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

    for seq in inputs {
        // First find the anchors, aka exact matches between
        // seqs and kmers in the index
        let mut seq_anchors: Vec<Anchor> = anchors_for_query(index, seq);
        
        // Chain close anchors together to find longer matches
        let seq_chains : Vec<Chain> = chain_anchors(&mut seq_anchors, index.kmer_length,bandwidth,
                                         max_gap, chain_min_n_anchors, secondary_chain_threshold,
                                         max_mismatch_rate, max_mapq);

        // Add chains to the Vec containing the previous chains
        chains.extend(seq_chains.into_iter());

        // Convert chains to strings and do the same
        let curr_chains_gaf_strings : Vec<String> = chains
            .iter()
            .map(|chain| write_chain_gaf(chain, index, &seq.name, seq.seq.len()))
            .collect();

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