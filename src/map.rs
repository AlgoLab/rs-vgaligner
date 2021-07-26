use crate::io::InputSequence;
use crate::index::Index;
use substring::Substring;

struct Anchor {
    query_begin : u64,
    query_end : u64,
    target_begin : u64,
    target_end : u64,
    max_chain_score : f64,
    best_predecessor : Option<Box<Anchor>>
}

impl Anchor {
    pub fn new(query_begin : u64, query_end : u64, target_begin : u64, target_end : u64) -> Self {
        Anchor {
            query_begin,
            query_end,
            target_begin,
            target_end,
            max_chain_score : 0f64,
            best_predecessor : None,
        }
    }
}

fn anchors_for_query(index : &Index, query : &InputSequence) -> Vec<Anchor> {
    let mut anchors : Vec<Anchor> = Vec::new();

    // Get kmers from the query, having the same size as the ones
    // stored in the index
    let kmer_size = index.kmer_length as usize;
    let query_kmers : Vec<String> = query.split_into_kmers(kmer_size);

    // Find anchors for each kmer in query
    for i in 0..query_kmers.len() {
        let kmer = query_kmers.get(i).unwrap();
        let ref_pos_vec = index.find_positions_for_query_kmer(&kmer.as_str());

        let mut anchors_for_kmer : Vec<Anchor> = Vec::new();
        for pos in ref_pos_vec {
            anchors_for_kmer.push(
                Anchor::new(i as u64, (i+kmer_size) as u64,
                pos.start, pos.end)
            );
        }

        anchors.extend(anchors_for_kmer);
    }

    anchors
}

// TODO: add other params
pub fn map_reads(index : &Index, inputs : &Vec<InputSequence>) {

    for seq in inputs {
        // First find the anchors, aka exact matches between
        // seqs and kmers in the index
        let anchors : Vec<Anchor> = anchors_for_query(index, seq);
    }

}