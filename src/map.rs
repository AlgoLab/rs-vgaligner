use crate::io::InputSequence;
use crate::index::Index;
use substring::Substring;
use std::cmp::min;
use itertools::Unique;
use std::collections::VecDeque;
use std::ops::Deref;
use bio::data_structures::interval_tree::*;
use core::cmp;
use std::cmp::Ordering::Equal;
use float_cmp::approx_eq;
use rayon::prelude::*;
use std::iter::FromIterator;

#[derive(Clone, Debug)]
pub struct Anchor {
    pub query_begin : u64,
    pub query_end : u64,
    target_begin : u64,
    target_begin_orient : bool,
    target_end : u64,
    target_end_orient : bool,
    max_chain_score : f64,
    best_predecessor : Option<Box<Anchor>>
}

impl Anchor {
    pub fn new(query_begin : u64, query_end : u64,
               target_begin : u64, target_begin_orient : bool,
               target_end : u64, target_end_orient : bool) -> Self {
        Anchor {
            query_begin,
            query_end,
            target_begin,
            target_begin_orient,
            target_end,
            target_end_orient,
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

        //println!("Curr kmer: {}", kmer);
        //println!("Positions on ref: {:#?}", ref_pos_vec);

        let mut anchors_for_kmer : Vec<Anchor> = Vec::new();
        for pos in ref_pos_vec {
            anchors_for_kmer.push(
                Anchor::new(i as u64, (i+kmer_size) as u64,
                pos.start, pos.start_orient,
                pos.end, pos.end_orient)
            );
        }

        anchors.extend(anchors_for_kmer);
    }

    anchors
}

#[derive(Debug, Clone)]
pub struct Chain {
    pub anchors : VecDeque<Anchor>,
    score : f64,
    mapping_quality : f64,
    is_secondary : bool,
    target_begin : u64,
    target_end : u64
}

impl PartialEq for Chain {
    fn eq(&self, other: &Self) -> bool {
        self.target_begin.eq(&other.target_begin)
            && self.target_end.eq(&other.target_end)
            && self.is_secondary.eq(&other.is_secondary)
            && approx_eq!(f64, self.score, other.score, ulps=2)
            && approx_eq!(f64, self.mapping_quality, other.mapping_quality, ulps=2)
    }
}

impl Eq for Chain {}

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

    pub fn compute_boundaries(&mut self, seed_length : u64, mismatch_rate : f64) {
        let first_target_end = self.anchors.front().unwrap().target_end;
        let first_target_end_orient = self.anchors.front().unwrap().target_end_orient;
        let last_target_begin = self.anchors.back().unwrap().target_begin;
        let last_target_begin_orient = self.anchors.back().unwrap().target_begin_orient;

        let first_target_begin = self.anchors.front().unwrap().target_begin;
        let first_target_begin_orient = self.anchors.front().unwrap().target_begin_orient;
        let last_target_end = self.anchors.back().unwrap().target_end;
        let last_target_end_orient = self.anchors.back().unwrap().target_end_orient;

        if first_target_begin_orient == last_target_end_orient &&
            first_target_begin < last_target_end &&
            self.score * (1f64+mismatch_rate) > (last_target_end - first_target_begin) as f64 {
            self.target_begin = first_target_begin;
            self.target_end = last_target_end;
        } else if first_target_end_orient == last_target_begin_orient &&
            first_target_end < last_target_begin {
            self.target_begin = first_target_end;
            self.target_end = last_target_begin;
        } else {
            self.score = - f64::MAX;
        }

    }
}


pub fn score_anchor(a : &Anchor, b : &Anchor, seed_length : &u64, max_gap : &u64) -> f64 {
    let mut score : f64 = a.max_chain_score;

    if a.query_end >= b.query_end ||
        !(((a.target_end_orient == b.target_end_orient) == a.target_begin_orient) == b.target_begin_orient) {
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
                false => 0.01f64 * (*seed_length as f64) * gap_length as f64 + 0.5f64 * f64::log2(gap_length as f64)
            };
            let match_length = min(min(query_length, target_length), *seed_length);
            score = f64::round((a.max_chain_score + match_length as f64 - gap_cost as f64) * 1000.0f64) / 1000.0f64 + query_overlap as f64;

        }
    }

    score
}

pub fn chains(anchors : &mut Vec<Anchor>, seed_length : u64, bandwidth : u64,
              max_gap : u64, chain_min_n_anchors : u64, secondary_chain_threshold : f64,
              mismatch_rate : f64, max_mapq : f64) -> Vec<Chain> {

    // First sort the anchors by their ending position
    anchors.sort_by(|a,b| a.target_end.cmp(&b.target_end));

    for i in 0..anchors.len() {

        let mut max_score = seed_length as f64;
        let mut best_predecessor : Option<Box<Anchor>> = None;

        let mut min_j : usize = 0;
        if bandwidth > i as u64 {
            min_j = 0;
        } else {
            min_j = i - bandwidth as usize;
        }

        for j in (i-1..min_j).rev() {
            let anchor_i = anchors.get(i).unwrap();
            let anchor_j = anchors.get(j).unwrap();

            let proposed_score = score_anchor(anchor_j, anchor_i, &seed_length, &max_gap);

            if proposed_score > max_score {
                max_score = proposed_score;
                best_predecessor = Some(Box::new(anchor_j.clone()));
            }
        }

        let anchor_i_mut = anchors.get_mut(i).unwrap();
        anchor_i_mut.max_chain_score = max_score;
        anchor_i_mut.best_predecessor = best_predecessor;

    }

    let mut chains: Vec<Chain> = Vec::new();


    // TODO: review this
    let mut i = anchors.len() - 1;
    while i >= 0 {
        let mut a = anchors.get_mut(i).unwrap();

        if a.best_predecessor.is_some() && a.max_chain_score > seed_length as f64 {
            let mut curr_chain = Chain::new();
            curr_chain.score = a.max_chain_score;

            loop {
                // First push the current anchor
                curr_chain.anchors.push_back(a.clone());

                // Then move backwards (predecessors)
                match a.best_predecessor {
                    None => break,
                    Some(ref mut b) => {
                        a = b;
                    }

                }
            }

            for anchor in &mut curr_chain.anchors {
                anchor.best_predecessor = None;
            }

            if curr_chain.anchors.len() >= chain_min_n_anchors as usize {
                curr_chain.anchors = VecDeque::from_iter(curr_chain.anchors.into_iter().rev());
                chains.push(curr_chain);
            }
        }
        i -= 1;
    }


    // Sort the chains by score in reverse order (higher first)
    chains.sort_by(|a,b| b.score.partial_cmp(&a.score).unwrap_or(Equal));

    // Create the Interval Tree
    let mut interval_tree : ArrayBackedIntervalTree<u64, Chain> = ArrayBackedIntervalTree::new();
    for chain in &chains {
        let query_begin = chain.anchors.front().unwrap().query_begin;
        let query_end = chain.anchors.back().unwrap().query_end;
        interval_tree.insert((query_begin..query_end), chain.clone());
    }
    interval_tree.index();

    // Find overlaps
    for mut chain in &mut chains {
        if !chain.processed() {
            let chain_begin = chain.anchors.front().unwrap().query_begin;
            let chain_end = chain.anchors.back().unwrap().query_end;
            let chain_length = chain_end - chain_begin;
            let mut best_secondary: Option<Chain> = None;

            let ovlp = interval_tree.find((chain_begin..chain_end));
            for value in &ovlp {
                let mut other_chain = value.data().clone();
                if other_chain != *chain && other_chain.score <= chain.score {
                    let other_begin = value.interval().start;
                    let other_end = value.interval().end;
                    let other_length = other_end - other_begin;
                    let ovlp_begin = cmp::max(chain_begin, other_begin);
                    let ovlp_end = cmp::min(chain_end, other_end);
                    let ovlp_length = ovlp_end - ovlp_begin;

                    if ovlp_length > other_length * secondary_chain_threshold as u64 {
                        other_chain.mapping_quality = 0f64;
                        other_chain.is_secondary = true;
                    }

                    if best_secondary.is_none() ||
                        best_secondary.is_some() && best_secondary.as_ref().unwrap().score < other_chain.score  {
                        best_secondary = Some(other_chain.clone());
                    }
                }
            }
            if best_secondary == None {
                chain.mapping_quality = max_mapq;
            } else {
                chain.mapping_quality =
                    40f64 * (1f64 -best_secondary.unwrap().score / chain.score)
                        * cmp::min(1, (chain.anchors.len()/10)) as f64
                        * f64::ln(chain.score);
            }
        }
    }


    for chain in &mut chains {
        chain.compute_boundaries(seed_length, mismatch_rate);
    }

    chains
}

/*
    GAF format
    https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf

    Col     Type    Description
    1       string  Query sequence name
    2       int     Query sequence length
    3       int     Query start (0-based; closed)
    4       int     Query end (0-based; open)
    5       char    Strand relative to the path: "+" or "-"

    6       string  Path matching /([><][^\s><]+(:\d+-\d+)?)+|([^\s><]+)/
    7       int     Path length
    8       int     Start position on the path (0-based)
    9       int     End position on the path (0-based)
    10      int     Number of residue matches
    11      int     Alignment block length
    12      int     Mapping quality (0-255; 255 for missing)
*/
pub fn write_chain_gaf(chain : &Chain, index : &Index,
                       query_name : &String, query_length : usize) -> String {

    let query_begin = chain.anchors.front().unwrap().query_begin;
    let query_end = chain.anchors.back().unwrap().query_end;
    let first_half_gaf_line = format!("{}\t{}\t{}\t{}\t{}\t",
                                      query_name,
                                      query_length,
                                      query_begin,
                                      query_end,
                                      "+");
    let path_length: u64 = 0;
    // TODO: continue this + needs fixes (original impl)
    let second_half_gaf_line = format!("\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                                       path_length,
                                       0,
                                       path_length,
                                       path_length,
                                       path_length,
                                       cmp::min(chain.mapping_quality as u64,254 as u64),
                                       "ta:Z:chain");

    format!("{}{}", first_half_gaf_line, second_half_gaf_line)
}

pub fn map_reads(index : &Index, inputs : &Vec<InputSequence>,
                 bandwidth : u64, max_gap : u64,
                 chain_min_n_anchors : u64, secondary_chain_threshold : f64,
                 max_mismatch_rate : f64, max_mapq : f64, write_chains : bool) {

    for seq in inputs {
        // First find the anchors, aka exact matches between
        // seqs and kmers in the index
        let mut anchors: Vec<Anchor> = anchors_for_query(index, seq);
        
        // Chain close anchors together to find longer matches
        let chains : Vec<Chain> = chains(&mut anchors, index.kmer_length,bandwidth,
                                         max_gap, chain_min_n_anchors, secondary_chain_threshold,
                                         max_mismatch_rate, max_mapq);

        if write_chains {
            for chain in &chains {
                write_chain_gaf(chain, index, &seq.name, seq.seq.len());
            }
        }
    }

}

#[cfg(test)]
mod test {
    use handlegraph::hashgraph::HashGraph;
    use handlegraph::mutablehandlegraph::MutableHandleGraph;
    use handlegraph::handle::Edge;
    use handlegraph::pathgraph::PathHandleGraph;
    use crate::index::Index;
    use crate::io::InputSequence;
    use crate::map::anchors_for_query;

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

    #[test]
    fn test_simple_anchors() {
        let mut graph = HashGraph::new();
        graph.create_handle("ACT".as_bytes(), 1);
        let mut index = Index::build(&graph, 3, 100,
                                     100, 7.0, None);
        let query : String = String::from("ACT");
        let input_seq = InputSequence::from_string(&query);
        let anchors = anchors_for_query(&index, &input_seq);
        assert_eq!(anchors.len(), 1);
    }

    #[test]
    fn test_anchors() {
        let mut graph = create_simple_graph();
        let index = Index::build(&graph, 3, 100,
                                     100, 7.0, None);
        let input_seq = InputSequence::from_string(&String::from("ACTGCA"));
        let anchors = anchors_for_query(&index, &input_seq);

        // There are at least 4 3-mers in the query, there can be more because
        // a kmer can appear in multiple positions
        assert!(anchors.len() >= 4)
    }

    #[test]
    fn test_no_anchors() {
        let mut graph = create_simple_graph();
        let index = Index::build(&graph, 3, 100,
                                 100, 7.0, None);
        let input_seq = InputSequence::from_string(&String::from("AAATTT"));
        let anchors = anchors_for_query(&index, &input_seq);
        assert_eq!(anchors.len(), 0)
    }


}