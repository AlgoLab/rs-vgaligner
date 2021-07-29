use crate::io::InputSequence;
use crate::index::Index;
use substring::Substring;
use std::cmp::min;
use itertools::Unique;
use std::collections::VecDeque;
use std::ops::Deref;

#[derive(Clone)]
struct Anchor {
    pub query_begin : u64,
    pub query_end : u64,
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

struct Chain {
    anchors : VecDeque<Anchor>,
    score : f64,
    mapping_quality : f64,
    is_secondary : bool,
    target_begin : u64,
    target_end : u64
}

impl Chain {
    pub fn new() -> Self {
        Chain {
            anchors : VecDeque::new(),
            score : 0f64,
            mapping_quality : f64::MIN,
            is_secondary : false,
            target_begin : 0,
            target_end : 0,
        }
    }

    pub fn processed(&self) -> bool {
        self.mapping_quality != f64::MIN
    }

    pub fn query_begin(&self) -> u64 {
        assert!(self.anchors.len() > 0);
        self.anchors.front().unwrap().query_begin
    }

    pub fn query_end(&self) -> u64 {
        assert!(self.anchors.len() > 0);
        self.anchors.back().unwrap().query_end
    }
}


pub fn score_anchor(a : &Anchor, b : &Anchor, seed_length : &u64, max_gap : &u64) -> f64 {
    let mut score : f64 = a.max_chain_score;
    // TODO: understand and add what's after ||
    if a.query_end >= b.query_end {
        score = -f64::MAX;
    } else {
        let query_length : u64 = min(b.query_begin-a.query_begin, b.query_end-a.query_end);

        let query_overlap = match a.query_end > b.query_end {
            true => a.query_end - b.query_end,
            false => 0
        };

        let target_length = min(b.target_begin-a.target_begin, b.target_end-a.target_end);

        // Abs of two unsigned integers doesn't make sense (in Rust)
        //let gap_length = (query_length - target_length).abs();
        // This should be equivalent
        let gap_length = match query_length > target_length {
            true => query_length - target_length,
            false => target_length - query_length
        };

        if gap_length > *max_gap {
            score = -f64::MAX;
        } else {
            let gap_cost = match gap_length == 0 {
                true => 0f64,
                // TODO: understand what this is
                false => 0.01f64 * (*seed_length as f64) * gap_length as f64 + 0.5f64 * f64::log2(gap_length as f64)
            };
            let match_length = min(min(query_length, target_length), *seed_length);
            // TODO: understand what this is
            score = f64::round((a.max_chain_score + match_length as f64 - gap_cost as f64) * 1000.0f64) / 1000.0f64 + query_overlap as f64;

        }
    }

    score
}

pub fn chains(anchors : &mut Vec<Anchor>, kmer_length : u64, seed_length : u64, bandwidth : u64, max_gap : u64) -> Vec<Chain> {

    // First sort the anchors by their ending position
    anchors.sort_by(|a,b| a.target_end.cmp(&b.target_end));

    for i in 0..anchors.len() {
        let anchor_i = anchors.get_mut(i).unwrap();
        anchor_i.max_chain_score = seed_length as f64;

        let mut min_j : usize = 0;
        if bandwidth > i as u64 {
            min_j = 0;
        } else {
            min_j = i - bandwidth as usize;
        }

        for j in (i-1..min_j).rev() {
            let anchor_j = anchors.get(j).unwrap();
            let proposed_score = score_anchor(&anchor_j, anchor_i, &seed_length, &max_gap);
            
            if proposed_score > anchor_i.max_chain_score {
                anchor_i.max_chain_score = proposed_score;
                // TODO: review this
                anchor_i.best_predecessor = Some(Box::new(anchor_j.clone()));
            }
        }

    }

    let chains : Vec<Chain> = Vec::new();
    let i = anchors.len() - 1;
    while i >= 0 {
        let a = anchors.get(i).unwrap();

        if a.best_predecessor.is_none() && a.max_chain_score > seed_length as f64 {
            let mut curr_chain = Chain::new();
            curr_chain.anchors.push_back(a.clone());
        }
    }



    chains
}

// TODO: add other params
pub fn map_reads(index : &Index, inputs : &Vec<InputSequence>, seed_length : u64, bandwidth : u64, max_gap : u64) {

    for seq in inputs {
        // First find the anchors, aka exact matches between
        // seqs and kmers in the index
        let mut anchors: Vec<Anchor> = anchors_for_query(index, seq);
        let chains : Vec<Chain> = chains(&mut anchors, index.kmer_length, seed_length, bandwidth, max_gap);
    }

}

#[cfg(test)]
mod test {
    use handlegraph::hashgraph::HashGraph;
    use handlegraph::mutablehandlegraph::MutableHandleGraph;
    use handlegraph::handle::Edge;
    use handlegraph::pathgraph::PathHandleGraph;

    /// This function creates a simple graph, used for debugging
    ///        | 2: CT \
    /// 1:  A            4: GCA
    ///        \ 3: GA |
    fn create_simple_graph() -> HashGraph {
        let mut graph: HashGraph = HashGraph::new();

        let h1 = graph.create_handle("A".as_bytes(), 1);
        let h2 = graph.create_handle("CT".as_bytes(), 2);
        let h3 = graph.create_handle("GA".as_bytes(), 3);
        let h4 = graph.create_handle("GCA".as_bytes(), 4);

        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h1, h3));
        graph.create_edge(&Edge(h2, h4));
        graph.create_edge(&Edge(h3, h4));

        let p1 = graph.create_path_handle("P1".as_bytes(), false);
        graph.append_step(&p1, h1);
        graph.append_step(&p1, h2);
        graph.append_step(&p1, h4);

        let p2 = graph.create_path_handle("P2".as_bytes(), false);
        graph.append_step(&p2, h1);
        graph.append_step(&p2, h3);
        graph.append_step(&p2, h4);

        graph
    }


}