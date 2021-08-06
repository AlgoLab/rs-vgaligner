use crate::index::Index;
use crate::io::InputSequence;
use crate::chain::{Anchor, Chain, write_chain_gaf, anchors_for_query, chain_anchors};

pub fn map_reads(index : &Index, inputs : &Vec<InputSequence>,
                 bandwidth : u64, max_gap : u64,
                 chain_min_n_anchors : u64, secondary_chain_threshold : f64,
                 max_mismatch_rate : f64, max_mapq : f64, write_chains : bool) {

    for seq in inputs {
        // First find the anchors, aka exact matches between
        // seqs and kmers in the index
        let mut anchors: Vec<Anchor> = anchors_for_query(index, seq);
        
        // Chain close anchors together to find longer matches
        let chains : Vec<Chain> = chain_anchors(&mut anchors, index.kmer_length,bandwidth,
                                         max_gap, chain_min_n_anchors, secondary_chain_threshold,
                                         max_mismatch_rate, max_mapq);

        if write_chains {
            for chain in &chains {
                write_chain_gaf(chain, index, &seq.name, seq.seq.len());
            }
        }
    }

}