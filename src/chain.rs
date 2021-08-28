use crate::io::InputSequence;
use crate::index::Index;
use substring::Substring;
use std::cmp::min;
use itertools::Unique;
use std::collections::VecDeque;
use std::ops::{Deref, Range};
use bio::data_structures::interval_tree::*;
use core::cmp;
use std::cmp::Ordering::Equal;
use float_cmp::approx_eq;
use rayon::prelude::*;
use std::iter::FromIterator;

// Instead of using pointers like in the original implementation,
// now each Anchor has an id, and the best predecessor is
// identified by its id
type anchor_id = u64;

/// An exact match between the input sequence and the graph. With respect to the Minimap2 paper,
/// where an anchor is a triple (x,y,w), in this implementation we have that:
/// - x = target_end
/// - y = query_end
/// - w = index.kmer_size
#[derive(Clone, Debug)]
pub struct Anchor {
    // The id of the anchor
    id : anchor_id,

    // These are positions on the query (i.e. seqs in the .fa/.fq file)
    pub query_begin : u64,
    pub query_end : u64,

    // These are positions on the graph linearization (position + fwd/rev)
    target_begin : u64,
    target_begin_orient : bool,
    target_end : u64,
    target_end_orient : bool,

    // Values used during chaining
    max_chain_score : f64,
    best_predecessor_id : Option<anchor_id>
}

impl Anchor {
    pub fn new(query_begin : u64, query_end : u64,
               target_begin : u64, target_begin_orient : bool,
               target_end : u64, target_end_orient : bool,
               id : u64) -> Self {
        Anchor {
            id,
            query_begin,
            query_end,
            target_begin,
            target_begin_orient,
            target_end,
            target_end_orient,
            max_chain_score : 0f64,
            best_predecessor_id : None,
        }
    }
}

/// Obtain all the anchors between the [index] and the [query] sequence
// Since the kmer positions are already stored in the index, one simply has to
// split the query into kmers, and then get the target positions from the index itself.
pub(crate) fn anchors_for_query(index : &Index, query : &InputSequence) -> Vec<Anchor> {
    let mut anchors : Vec<Anchor> = Vec::new();

    // Get kmers from the query, having the same size as the ones
    // stored in the index
    let kmer_size = index.kmer_length as usize;
    let query_kmers : Vec<String> = query.split_into_kmers(kmer_size);

    let mut id : anchor_id = 0;
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
                            pos.end, pos.end_orient, id)
            );
            id += 1;
        }

        anchors.extend(anchors_for_kmer);
    }

    anchors
}

/// Multiple anchors chained together to get a longer match
#[derive(Debug, Clone)]
pub struct Chain {
    // A Chain is made up by one or more anchors
    pub anchors : VecDeque<Anchor>,
    score : f64,
    mapping_quality : f64,
    is_secondary : bool,
    target_begin : u64,
    target_begin_orient : bool,
    target_end : u64,
    target_end_orient : bool
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
            target_begin_orient: true,
            target_end : 0,
            target_end_orient : true
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

    // find the longest contiguous range covered by the chain
    // where the size is no more than some function of our seed_length * number of anchors * 1+mismatch_rate
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

        // TODO: review this
        // Removed the following line because it would overflow (i.e. a-b when a < b)
        //let target_length = min(b.target_begin-a.target_begin, b.target_end-a.target_end);

        let target_begin_diff = match b.target_begin > a.target_begin {
            true => b.target_begin - a.target_begin,
            false => a.target_begin - b.target_begin
        };
        let target_end_diff = match b.target_end > a.target_end {
            true => b.target_end - a.target_end,
            false => a.target_end - b.target_end
        };
        let target_length = min(target_begin_diff, target_end_diff);

        // First we compute the inner part of gamma_c, aka the "gap length".
        // RUST DETAIL: Abs of two unsigned integers doesn't make sense (in Rust)
        //let gap_length = (query_length - target_length).abs();
        // The following match should be equivalent.
        let gap_length = match query_length > target_length {
            true => query_length - target_length,
            false => target_length - query_length
        };

        if gap_length > *max_gap {
            score = -f64::MAX;
        } else {
            // This is gamma_c(l) in the paper, aka the "gap cost". This is also equivalent
            // to B(j,i).
            let gap_cost = match gap_length == 0 {
                true => 0f64,
                false => 0.01f64 * (*seed_length as f64) * gap_length as f64 + 0.5f64 * f64::log2(gap_length as f64)
            };

            // This is a(j,i)
            let match_length = min(min(query_length, target_length), *seed_length);

            // This is f(j) + a(j,i) - B(j,i)
            score = f64::round((a.max_chain_score + match_length as f64 - gap_cost as f64) * 1000.0f64) / 1000.0f64 + query_overlap as f64;
        }
    }

    score
}

pub fn chain_anchors(anchors : &mut Vec<Anchor>, seed_length : u64, bandwidth : u64,
              max_gap : u64, chain_min_n_anchors : u64, secondary_chain_threshold : f64,
              mismatch_rate : f64, max_mapq : f64) -> Vec<Chain> {

    // ----- STEP 1 : finding the optimal chaining scores -----

    // First sort the anchors by their ending position
    anchors.par_sort_by(|a,b| a.target_end.cmp(&b.target_end));

    // Then, compute the maximal chaining score up to the current anchor.
    // This score is called f(i) in the paper.
    // NOTE: Starts at 1 because anchors[0] has no previous values
    for i in 1..anchors.len() {

        let mut min_j = match bandwidth > i as u64 {
            true => 0,
            false => i - bandwidth as usize,
        };

        //RUST NOTE: in a range (a..b), a MUST ALWAYS be less than b (even when going backwards
        // with .rev()!) otherwise it will be empty. Previously, due to this error,
        // the following loop would not work.

        // RUST NOTE 2: the usage of references here a bit of a hack...
        // I should get a mutable reference to anchors[i] in this scope
        // and a non mutable one to anchors[j] inside the nested loop,
        // but Rust does not like that. TODO (again)?
        for j in (min_j..i-1).rev() {
            let anchor_j = anchors.get(j).unwrap().clone();
            let anchor_i = anchors.get_mut(i).unwrap();

            // This is where we compute f(j) + a(j,i) - B(j,i)
            let proposed_score = score_anchor(&anchor_j, anchor_i, &seed_length, &max_gap);

            // We are now obtaining f(i)
            if proposed_score > anchor_i.max_chain_score {
                anchor_i.max_chain_score = proposed_score;
                anchor_i.best_predecessor_id = Some(anchor_j.id);
            }

        }

    }


    // ----- STEP 2: Finding all the chains with no anchors used in multiple chains (backtracking) -----
    let mut chains: Vec<Chain> = Vec::new();

    if !anchors.is_empty() {
        let mut i = anchors.len() - 1;

        loop {
            let mut a = anchors.get_mut(i).unwrap();

            if a.best_predecessor_id.is_some() && a.max_chain_score > seed_length as f64 {

                // Mark predecessor as None (-> a won't be used again)
                let mut pred_id = a.best_predecessor_id.clone();
                a.best_predecessor_id = None;

                let mut curr_chain = Chain::new();
                curr_chain.anchors.push_back(a.clone());
                curr_chain.score = a.max_chain_score;

                // Go back the best predecessors until a None is found
                loop {
                    match pred_id {
                        Some(b) => {
                            // Get the anchor relative to the best predecessor
                            let mut pred_anchor = anchors.get_mut(b as usize).unwrap();

                            // Check if it has already been used (= its best predecessor is None)
                            if pred_anchor.best_predecessor_id.is_none() {
                                pred_id = None;
                            } else {
                                pred_id = pred_anchor.best_predecessor_id.clone();
                                pred_anchor.best_predecessor_id = None;
                                curr_chain.anchors.push_back(pred_anchor.clone());
                            }
                        }
                        None => break,
                    }
                }

                if curr_chain.anchors.len() >= chain_min_n_anchors as usize {
                    curr_chain.anchors = VecDeque::from_iter(curr_chain.anchors.into_iter().rev());
                    chains.push(curr_chain.clone());
                }
            }

            // Decrease counter for next iteration
            if i>0 {
                i -= 1;
            } else {
                break;
            }

        }
    }

    // ----- STEP 3: Identifying primary chains (= chains with little or no overlap on the query) -----

    // Sort the chains by score in reverse order (higher first)
    chains.par_sort_by(|a,b| b.score.partial_cmp(&a.score).unwrap_or(Equal));

    // Create the Interval Tree
    let mut interval_tree : ArrayBackedIntervalTree<u64, Chain> = ArrayBackedIntervalTree::new();
    for chain in &chains {
        let query_begin = chain.anchors.front().unwrap().query_begin;
        let query_end = chain.anchors.back().unwrap().query_end;

        // TODO: review this, maybe query_begin should always be less than query_end?
        if query_begin < query_end {
            interval_tree.insert((query_begin..query_end), chain.clone());
        }

    }
    interval_tree.index();

    // Find overlaps
    for mut chain in &mut chains {
        // TODO: same as line 333
        if !chain.processed() && chain.query_begin() < chain.query_end() {
            let chain_begin = chain.anchors.front().unwrap().query_begin;
            let chain_end = chain.anchors.back().unwrap().query_end;
            //let chain_length = chain_end - chain_begin;
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
pub fn extract_graph_po_range(index : &Index, chains : &Vec<Chain>, query : &InputSequence) -> Vec<Range<u64>> {
    let mut start_end_vec : Vec<Range<u64>> = Vec::new();

    for chain in chains {
        let min_node = chain.anchors.par_iter()
                                    .map(|a| index.get_node_from_pos(a.target_begin as usize, a.target_begin_orient))
                                    .min().unwrap();

        let max_node = chain.anchors.par_iter()
                                    .map(|a| index.get_node_from_pos(a.target_end as usize, a.target_end_orient))
                                    .max().unwrap();

        start_end_vec.push(min_node..max_node);
    }

    start_end_vec
}
 */

/// Find the leftmost starting node(handle) and rightmost node(handle) for each chain
pub fn extract_graph_po_range(index : &Index, chains : &Vec<Chain>) -> Vec<Range<u64>> {
    let mut start_end_vec : Vec<Range<u64>> = chains.par_iter()
                                                    .map(|c| find_range_single_chain(index, c))
                                                    .collect();

    start_end_vec
}

fn find_range_single_chain(index : &Index, chain : &Chain) -> Range<u64> {
    let min_node = chain.anchors.par_iter()
        .map(|a| index.get_node_from_pos(a.target_begin as usize, a.target_begin_orient))
        .min().unwrap();

    let max_node = chain.anchors.par_iter()
        .map(|a| index.get_node_from_pos(a.target_end as usize, a.target_end_orient))
        .max().unwrap();

    (min_node..max_node)
}

pub fn extract_query_subsequence(index : &Index, chains : &Vec<Chain>, query : &InputSequence) -> Vec<String> {
    let subseqs : Vec<String> = chains.par_iter()
                                      .map(|c| extract_single_subsequence(index,
                                                                          c.target_begin as usize, c.target_begin_orient,
                                                                          c.target_end as usize, c.target_end_orient))
                                      .collect();
    subseqs
}

fn extract_single_subsequence(index: &Index, begin: usize, begin_orient: bool, end: usize, end_orient: bool) -> String {
    let substring = match (begin_orient, end_orient) {
        (true, true) => index.seq_fwd.substring(begin, end).to_string(),
        (false, false) => index.seq_rev.substring(begin, end).to_string(),

        // TODO: this is 100% not right, maybe I should take a part from fwd and another from rev?
        _ => index.seq_fwd.substring(begin, end).to_string(),

    };

    substring
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
/// Convert a chain to a string in GAF (Graph Alignment Format)
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
                                       "",
                                       path_length,
                                       path_length,
                                       path_length,
                                       path_length,
                                       cmp::min(chain.mapping_quality as u64,254 as u64),
                                       "ta:Z:chain");

    format!("{}{}", first_half_gaf_line, second_half_gaf_line)
}

#[cfg(test)]
mod test {
    use handlegraph::hashgraph::HashGraph;
    use handlegraph::mutablehandlegraph::MutableHandleGraph;
    use handlegraph::handle::Edge;
    use handlegraph::pathgraph::PathHandleGraph;
    use crate::index::Index;
    use crate::io::InputSequence;
    use crate::chain::{anchors_for_query, chain_anchors, Chain};
    use gfa::parser::GFAParser;
    use gfa::gfa::GFA;
    use std::path::PathBuf;

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

    #[test]
    fn test_no_anchors_2() {
        let mut graph = create_simple_graph();
        let index = Index::build(&graph, 3, 100,
                                 100, 7.0, None);
        let input_seq = InputSequence::from_string(&String::from(""));
        let anchors = anchors_for_query(&index, &input_seq);
        assert_eq!(anchors.len(), 0)
    }

    #[test]
    fn test_chains() {
        let mut graph = create_simple_graph();
        let index = Index::build(&graph, 3, 100,
                                 100, 7.0, None);
        let input_seq = InputSequence::from_string(&String::from("ACTGCA"));
        let mut anchors = anchors_for_query(&index, &input_seq);

        let chains: Vec<Chain> = chain_anchors(&mut anchors, index.kmer_length, 50,
                                               1000, 1, 0.5f64,
                                               0.1f64, 60.0f64);
        assert!(chains.is_empty());
        //println!("Chains: {:#?}", chains);
        //println!("Chains length: {}", chains.len());
    }

    #[test]
    fn test_chains_2() {

        // Create HashGraph from GFA
        let parser = GFAParser::new();
        let gfa: GFA<usize, ()> = parser.parse_file(&PathBuf::from("./test/test.gfa")).unwrap();
        let graph = HashGraph::from_gfa(&gfa);

        let index = Index::build(&graph, 11, 100,
                                 100, 7.0, None);
        let input_seq = InputSequence::from_string(&String::from(index.seq_fwd.clone()));
        let mut anchors = anchors_for_query(&index, &input_seq);

        //println!("Anchors len: {}", anchors.len());
        //println!("Anchors: {:#?}", anchors);
        assert!(!anchors.is_empty());

        let chains: Vec<Chain> = chain_anchors(&mut anchors, index.kmer_length, 50,
                                               1000, 2, 0.5f64,
                                               0.1f64, 60.0f64);
        assert!(!chains.is_empty());
        //println!("Chains_2: {:#?}", chains);
        //println!("Chains_2 length: {}", chains.len());
    }

    #[test]
    fn test_no_chains() {
        let mut graph = create_simple_graph();
        let index = Index::build(&graph, 3, 100,
                                 100, 7.0, None);
        let input_seq = InputSequence::from_string(&String::from("AAATTT"));
        let mut anchors = anchors_for_query(&index, &input_seq);
        let chains: Vec<Chain> = chain_anchors(&mut anchors, index.kmer_length, 50,
                                               1000, 3, 0.5f64,
                                               0.1f64, 60.0f64);
        assert!(chains.is_empty());
    }

}