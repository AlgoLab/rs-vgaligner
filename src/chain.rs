use core::cmp;
use std::cmp::{min, Ordering};
use std::cmp::Ordering::Equal;
use std::collections::VecDeque;
use std::iter::FromIterator;
use std::ops::Range;

use bio::data_structures::interval_tree::*;
use float_cmp::approx_eq;
use serde::de::Unexpected::Seq;
//use rayon::prelude::*;

use crate::index::Index;
use crate::io::QuerySequence;
use crate::kmer::{SeqOrient, SeqPos};
use crate::kmer::SeqOrient::Forward;

// Instead of using pointers like in the original implementation,
// now each Anchor has an id, and the best predecessor is
// identified by its id
type AnchorId = u64;

/// An exact match between the input sequence and the graph. With respect to the Minimap2 paper,
/// where an anchor is a triple (x,y,w), in this implementation we have that:
/// - x = target_end
/// - y = query_end
/// - w = index.kmer_size
#[derive(Clone, Debug)]
pub struct Anchor {
    // The id of the anchor
    pub id: AnchorId,

    // These are positions on the query (i.e. seqs in the .fa/.fq file)
    pub query_begin: u64,
    pub query_end: u64,

    // These are positions on the graph linearization (position + fwd/rev)
    pub target_begin: SeqPos,
    pub target_end: SeqPos,

    // Values used during chaining
    pub max_chain_score: f64,
    pub best_predecessor_id: Option<AnchorId>,
}
impl Anchor {
    pub fn new(
        query_begin: u64,
        query_end: u64,
        target_begin: SeqPos,
        target_end: SeqPos,
        id: u64,
        max_chain_score: f64,
    ) -> Self {
        Anchor {
            id,
            query_begin,
            query_end,
            target_begin,
            target_end,
            max_chain_score,
            best_predecessor_id: None,
        }
    }

    pub fn get_end_seqpos_inclusive(&self) -> SeqPos {
        let mut end_seqpos = self.target_end.clone();
        // Remove 1 because rhs non-inclusive, we want the actual end position
        end_seqpos.position -= 1;
        end_seqpos
    }

    pub fn get_detailed_graph_position(&self, index: &Index) -> AnchorPosOnGraph  {
        AnchorPosOnGraph::new(self, index)
    }
}

pub struct AnchorPosOnGraph {
    // Anchor start
    pub start_orient: SeqOrient,
    pub start_node: u64,
    pub start_offset_from_node: u64,

    // Anchor end
    pub end_orient: SeqOrient,
    pub end_node: u64,
    pub end_offset_from_node: u64,
}

impl AnchorPosOnGraph {
    pub fn new(anchor: &Anchor, index: &Index) -> AnchorPosOnGraph {
        let first_handle = index.handle_from_seqpos(&anchor.target_begin);
        let first_node_start = index.get_bv_select(u64::from(first_handle.id()));
        let (first_node_orient, first_node_offset) = match first_handle.is_reverse() {
            false => (SeqOrient::Forward, anchor.target_begin.position - first_node_start as u64),
            //true => (SeqOrient::Reverse, first_node_start as u64 - anchor.target_begin.position)
            true => (SeqOrient::Reverse, anchor.target_begin.position - first_node_start as u64)
        };

        let last_handle = index.handle_from_seqpos(&anchor.get_end_seqpos_inclusive());
        let last_node_start = index.get_bv_select(u64::from(last_handle.id()));
        let (last_node_orient, last_node_offset) = match last_handle.is_reverse() {
            false => (SeqOrient::Forward, anchor.target_end.position - last_node_start as u64),
            //true => (SeqOrient::Reverse, last_node_start as u64 - anchor.target_end.position)
            true => (SeqOrient::Reverse, anchor.target_end.position - last_node_start as u64)
        };

        AnchorPosOnGraph {
            start_orient: first_node_orient,
            start_node: first_handle.unpack_number(),
            start_offset_from_node: first_node_offset,
            end_orient: last_node_orient,
            end_node: last_handle.unpack_number(),
            end_offset_from_node: last_node_offset
        }
    }
}

/// Obtain all the anchors between the [index] and the [query] sequence
// Since the kmer positions are already stored in the index, one simply has to
// split the query into kmers, and then get the target positions from the index itself.
// NOTE: if a kmer appears in multiple positions, ALL the resulting anchors will be emitted
pub(crate) fn anchors_for_query(index: &Index, query: &QuerySequence, only_forward: bool) -> Vec<Anchor> {
    let mut anchors: Vec<Anchor> = Vec::new();

    // Get kmers from the query, having the same size as the ones
    // stored in the index
    let kmer_size = index.kmer_length as usize;
    let query_kmers: Vec<String> = query.split_into_kmers(kmer_size);

    let mut id: AnchorId = 0;
    // Find anchors for each kmer in query
    for i in 0..query_kmers.len() {
        let kmer = query_kmers.get(i).unwrap();
        let ref_pos_vec = index.find_positions_for_query_kmer(&kmer.as_str());

        let mut anchors_for_kmer: Vec<Anchor> = Vec::new();
        for pos in ref_pos_vec {
            if (only_forward && pos.start.orient == Forward && pos.end.orient == Forward) || (!only_forward)  {
                anchors_for_kmer.push(Anchor::new(
                    i as u64,
                    (i + kmer_size) as u64,
                    pos.start,
                    pos.end,
                    id,
                    index.kmer_length as f64,
                ));
                id += 1;
            }
        }

        anchors.extend(anchors_for_kmer);
    }

    anchors
}

/// Multiple anchors chained together to get a longer match
#[derive(Debug, Clone)]
pub struct Chain {
    // A Chain is made up by one or more anchors
    pub anchors: VecDeque<Anchor>,
    pub score: f64,
    pub mapping_quality: f64,
    pub is_secondary: bool,
    pub target_begin: SeqPos,
    pub target_end: SeqPos,
    pub query: QuerySequence,
    pub is_placeholder: bool,
}
impl PartialEq for Chain {
    fn eq(&self, other: &Self) -> bool {
        self.target_begin.eq(&other.target_begin)
            && self.target_end.eq(&other.target_end)
            && self.is_secondary.eq(&other.is_secondary)
            && approx_eq!(f64, self.score, other.score, ulps = 2)
            && approx_eq!(f64, self.mapping_quality, other.mapping_quality, ulps = 2)
    }
}
impl Eq for Chain {}
impl Chain {
    pub fn new() -> Self {
        Chain {
            anchors: VecDeque::new(),
            score: 0f64,
            mapping_quality: f64::MIN,
            is_secondary: false,
            target_begin: SeqPos::new(SeqOrient::Forward, 0),
            target_end: SeqPos::new(SeqOrient::Forward, 0),
            query: QuerySequence::new(),
            is_placeholder: false,
        }
    }

    pub fn new_placeholder() -> Self {
        Chain {
            anchors: VecDeque::new(),
            score: 0f64,
            mapping_quality: f64::MIN,
            is_secondary: false,
            target_begin: SeqPos::new(SeqOrient::Forward, 0),
            target_end: SeqPos::new(SeqOrient::Forward, 0),
            query: QuerySequence::new(),
            is_placeholder: true,
        }
    }

    pub fn processed(&self) -> bool {
        self.mapping_quality != f64::MIN
    }

    pub fn query_begin(&self) -> u64 {
        assert!(!self.anchors.is_empty());
        self.anchors.front().unwrap().query_begin
    }

    pub fn query_end(&self) -> u64 {
        assert!(!self.anchors.is_empty());
        self.anchors.back().unwrap().query_end
    }

    // find the longest contiguous range covered by the chain
    // where the size is no more than some function of our seed_length * number of anchors * 1+mismatch_rate
    pub fn compute_boundaries(&mut self, _seed_length: u64, mismatch_rate: f64) {
        let first_target_end = self.anchors.front().unwrap().target_end;
        let last_target_begin = self.anchors.back().unwrap().target_begin;

        let first_target_begin = self.anchors.front().unwrap().target_begin;
        let last_target_end = self.anchors.back().unwrap().target_end;

        if first_target_begin.orient == last_target_end.orient
            && first_target_begin < last_target_end
            && self.score * (1f64 + mismatch_rate)
                > (last_target_end.position - first_target_begin.position) as f64
        {
            self.target_begin = first_target_begin.clone();
            self.target_end = last_target_end.clone();
        } else if first_target_end.orient == last_target_begin.orient
            && first_target_end < last_target_begin
        {
            self.target_begin = first_target_end.clone();
            self.target_end = last_target_begin.clone();
        } else {
            self.score = -f64::MAX;
        }
    }

    pub fn find_query_start_end(&self) -> Range<u64> {
        let min_query_pos: u64 = self.anchors.iter().map(|a| a.query_begin).min().unwrap();

        let max_query_pos: u64 = self.anchors.iter().map(|a| a.query_end).max().unwrap();

        min_query_pos..max_query_pos
    }
}

pub fn score_anchor(a: &Anchor, b: &Anchor, seed_length: &u64, max_gap: &u64) -> f64 {
    let mut score: f64 = a.max_chain_score;

    if a.query_end >= b.query_end
    //if a.query_end > b.query_end
        || (a.target_end.orient == b.target_end.orient && a.target_end.position >= b.target_end.position)
        || !(a.target_end.orient == b.target_end.orient
            && a.target_begin.orient == b.target_begin.orient
            && a.target_end.orient == b.target_begin.orient
            && a.target_begin.orient == b.target_end.orient)
    {
        println!("Min score between i: {:?} (orients: {:?}/{:?}) and j: {:?} (orients: {:?}/{:?})",
                 b.id, b.target_begin.orient, b.target_end.orient,
                 a.id, a.target_begin.orient, a.target_end.orient);

        if a.query_end >= b.query_end {
            println!("A ends after B in the query ({} vs {})", a.query_end, b.query_end);
        } else if a.target_end.orient == b.target_end.orient && a.target_end.position >= b.target_end.position {
            println!("A ends after B in the linearization ({} vs {})", a.target_end.position, b.target_end.position);
        } else {
            println!("A and B have different orients");
        }

        score = -f64::MAX;
    } else {
        let query_length: u64 = min(b.query_begin - a.query_begin, b.query_end - a.query_end);

        let query_overlap = match a.query_end > b.query_end {
            true => a.query_end - b.query_end,
            false => 0,
        };

        // TODO: review this
        // Removed the following line because it would overflow (i.e. a-b when a < b)
        //let target_length = min(b.target_begin-a.target_begin, b.target_end-a.target_end);

        let target_begin_diff = match b.target_begin > a.target_begin {
            true => b.target_begin.position - a.target_begin.position,
            false => a.target_begin.position - b.target_begin.position,
        };
        let target_end_diff = match b.target_end > a.target_end {
            true => b.target_end.position - a.target_end.position,
            false => a.target_end.position - b.target_end.position,
        };
        let target_length = min(target_begin_diff, target_end_diff);

        // First we compute the inner part of gamma_c, aka the "gap length".
        // RUST DETAIL: Abs of two unsigned integers doesn't make sense (in Rust)
        //let gap_length = (query_length - target_length).abs();
        // The following match should be equivalent.
        let gap_length = match query_length > target_length {
            true => query_length - target_length,
            false => target_length - query_length,
        };

        if gap_length > *max_gap {
            score = -f64::MAX;
        } else {
            // This is gamma_c(l) in the paper, aka the "gap cost". This is also equivalent
            // to B(j,i).
            let gap_cost = match gap_length == 0 {
                true => 0f64,
                false => {
                    0.01f64 * (*seed_length as f64) * gap_length as f64
                        + 0.5f64 * f64::log2(gap_length as f64)
                }
            };

            // This is a(j,i)
            let match_length = min(min(query_length, target_length), *seed_length);

            // This is f(j) + a(j,i) - B(j,i)
            score =
                f64::round((a.max_chain_score + match_length as f64 - gap_cost as f64) * 1000.0f64)
                    / 1000.0f64
                    + query_overlap as f64;
        }
    }

    score
}

pub fn chain_anchors(
    anchors: &mut Vec<Anchor>,
    seed_length: u64,
    bandwidth: u64,
    max_gap: u64,
    chain_min_n_anchors: u64,
    secondary_chain_threshold: f64,
    mismatch_rate: f64,
    max_mapq: f64,
    query: &QuerySequence,
) -> Vec<Chain> {
    // ----- STEP 1 : finding the optimal chaining scores -----

    // First sort the anchors by their ending position
    //anchors.sort_by(|a, b| a.target_end.position.cmp(&b.target_end.position));

    anchors.sort_by(|a,b|
        match b.target_end.orient.cmp(&a.target_end.orient) {
            Ordering::Equal => a.target_end.position.cmp(&b.target_end.position),
            other => other,
        });

    println!("Sorted anchors: {:?}", anchors.iter().map(|x| x.id).collect::<Vec<u64>>());
    //println!("Query: {}", query.seq);
    //println!("Query: {}, Anchors: {:#?}", query.seq, anchors);

    // Then, compute the maximal chaining score up to the current anchor.
    // This score is called f(i) in the paper.
    // NOTE: Starts at 1 because anchors[0] has no previous values
    for i in 1..anchors.len() {

        let min_j = match bandwidth > i as u64 {
            true => 0,
            false => i - bandwidth as usize,
        };

        // RUST NOTE: in a range (a..b), a MUST ALWAYS be less than b (even when going backwards
        // with .rev()!) otherwise it will be empty. Previously, due to this error,
        // the following loop would not work.

        // RUST NOTE 2: the usage of references here a bit of a hack...
        // I should get a mutable reference to anchors[i] in this scope
        // and a non mutable one to anchors[j] inside the nested loop,
        // but Rust does not like that. TODO (again)?
        for j in (min_j..i).rev() {
            let anchor_j = anchors.get(j).unwrap().clone();
            let anchor_i = anchors.get_mut(i).unwrap();

            // This is where we compute f(j) + a(j,i) - B(j,i)
            let proposed_score = score_anchor(&anchor_j, anchor_i, &seed_length, &max_gap);

            println!("For i {} and j {} proposed score is: {}", anchor_i.id, anchor_j.id, proposed_score);

            // We are now obtaining f(i)
            if proposed_score > anchor_i.max_chain_score {
                println!("Proposed: {} vs max_chain: {} (predecessor: {:?})", proposed_score, anchor_i.max_chain_score, anchor_i.best_predecessor_id);

                anchor_i.max_chain_score = proposed_score;
                anchor_i.best_predecessor_id = Some(anchor_j.id);

                println!("Updated anchor {}'s score to: {} (predecessor: {:?})\n", anchor_i.id, anchor_i.max_chain_score, anchor_i.best_predecessor_id);
            }
        }
        println!("\n");
    }

    // ----- STEP 2: Finding all the chains with no anchors used in multiple chains (backtracking) -----
    let mut chains: Vec<Chain> = Vec::new();

    if !anchors.is_empty() {

        for anchor in anchors.iter() {
            println!("For i {} predecessor is {:?}", anchor.id, anchor.best_predecessor_id);
        }
        println!("\n");

        let mut i = anchors.len() - 1;

        loop {

            let mut curr_anchor = anchors.get_mut(i).unwrap();

            if curr_anchor.best_predecessor_id.is_some() && curr_anchor.max_chain_score > seed_length as f64 {

                // Create the chain
                let mut curr_chain = Chain::new();
                curr_chain.query = query.clone();

                // Backtrack until an anchor with no predecessors is found
                while let Some(pred_id) = curr_anchor.best_predecessor_id {
                    // Set the current anchor's predecessor to None, to avoid re-using it in the future
                    curr_anchor.best_predecessor_id = None;

                    // Add the current anchor
                    curr_chain.anchors.push_back(curr_anchor.clone());
                    println!("Curr chain: {:?}", curr_chain.anchors.iter().map(|x| x.id).collect::<Vec<u64>>());

                    // Get the position of the previous anchor in the ordering
                    // NOTE: the anchor's id and the position in the anchors vec may be different
                    let pred_pos = anchors.iter().position(|x| x.id == pred_id).unwrap();

                    // Move the "pointer" to the predecessor
                    curr_anchor = anchors.get_mut(pred_pos).unwrap();
                }
                curr_chain.anchors.push_back(curr_anchor.clone());

            /*
            let mut a = anchors.get_mut(i).unwrap();
            if a.best_predecessor_id.is_some() && a.max_chain_score > seed_length as f64 {
                // Mark predecessor as None (-> a won't be used again)
                let mut pred_id = a.best_predecessor_id;
                a.best_predecessor_id = None;

                // Create chain + add the query it refers to
                let mut curr_chain = Chain::new();
                curr_chain.query = query.clone();

                curr_chain.anchors.push_back(a.clone());
                curr_chain.score = a.max_chain_score;

                // Go back the best predecessors until a None is found
                loop {
                    match pred_id {
                        Some(b) => {
                            // Get the anchor relative to the best predecessor
                            let mut pred_anchor = anchors.get_mut(b as usize).unwrap();
                            curr_chain.anchors.push_back(pred_anchor.clone());
                            pred_id = pred_anchor.best_predecessor_id;
                            pred_anchor.best_predecessor_id = None;
                            println!("Curr chain: {:?}", curr_chain.anchors.iter().map(|x| x.id).collect::<Vec<u64>>());

                            /*
                            // Check if it has already been used (= its best predecessor is None)
                            if pred_anchor.best_predecessor_id.is_none() {
                                curr_chain.anchors.push_back(pred_anchor.clone());
                                pred_id = None;
                            } else {
                                pred_id = pred_anchor.best_predecessor_id;
                                pred_anchor.best_predecessor_id = None;
                                curr_chain.anchors.push_back(pred_anchor.clone());
                                println!("Curr chain: {:?}", curr_chain.anchors.iter().map(|x| x.id).collect::<Vec<u64>>());
                                //println!("Pred id: {:?}", pred_id);
                            }
                             */
                        }
                        None => break,
                    }
                }

             */

                if curr_chain.anchors.len() >= chain_min_n_anchors as usize {
                    curr_chain.anchors = VecDeque::from_iter(curr_chain.anchors.into_iter().rev());
                    chains.push(curr_chain.clone());
                }
            }

            // Decrease counter for next iteration
            if i > 0 {
                i -= 1;
            } else {
                break;
            }
        }
    }

    // ----- STEP 3: Identifying primary chains (= chains with little or no overlap on the query) -----

    // Sort the chains by score in reverse order (higher first)
    chains.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap_or(Equal));
    /*
    chains.sort_by(|a,b| match b.anchors.len().cmp(&a.anchors.len()) {
        Ordering::Equal => b.score.partial_cmp(&a.score).unwrap_or(Equal),
        other => other
    });
     */

    //println!("Anchors per chain: {:#?}", chains.iter().map(|c| c.anchors.len()).collect::<Vec<usize>>());
    //println!("Chains are: {:#?}", chains);
    for (i, chain) in chains.iter().enumerate() {
        println!("Chain {} has anchors: {:?}", i, chain.anchors.iter().map(|x| x.id).collect::<Vec<u64>>());
    }

    /*
    // Create the Interval Tree
    let mut interval_tree: ArrayBackedIntervalTree<u64, Chain> = ArrayBackedIntervalTree::new();
    for chain in &chains {
        let query_begin = chain.anchors.front().unwrap().query_begin;
        let query_end = chain.anchors.back().unwrap().query_end;

        // TODO: review this, maybe query_begin should always be less than query_end?
        if query_begin < query_end {
            interval_tree.insert(query_begin..query_end, chain.clone());
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

            let ovlp = interval_tree.find(chain_begin..chain_end);
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

                    if best_secondary.is_none()
                        || best_secondary.is_some()
                            && best_secondary.as_ref().unwrap().score < other_chain.score
                    {
                        best_secondary = Some(other_chain.clone());
                    }
                }
            }
            if best_secondary == None {
                chain.mapping_quality = max_mapq;
            } else {
                chain.mapping_quality = 40f64
                    * (1f64 - best_secondary.unwrap().score / chain.score)
                    * cmp::min(1, chain.anchors.len() / 10) as f64
                    * f64::ln(chain.score);
            }
        }
    }

    for chain in &mut chains {
        chain.compute_boundaries(seed_length, mismatch_rate);
    }
     */

    if chains.is_empty() {
        // Create a placeholder chain -- i.e. no alignment was found
        let mut placeholder_chain = Chain::new_placeholder();
        placeholder_chain.query = query.clone();
        chains.push(placeholder_chain)
    }

    // There should always be at least one chain -- either a true one or a placeholder one
    assert!(!chains.is_empty());

    chains
}

/// Convert a chain to a string in GAF (Graph Alignment Format)
pub fn write_chain_gaf(
    chain: &Chain,
    _index: &Index,
    query_name: &String,
    query_length: usize,
) -> String {
    let query_begin = chain.anchors.front().unwrap().query_begin;
    let query_end = chain.anchors.back().unwrap().query_end;
    let first_half_gaf_line = format!(
        "{}\t{}\t{}\t{}\t{}\t",
        query_name, query_length, query_begin, query_end, "+"
    );
    let path_length: u64 = 0;
    // TODO: continue this + needs fixes (original impl)
    let second_half_gaf_line = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        "",
        path_length,
        path_length,
        path_length,
        path_length,
        cmp::min(chain.mapping_quality as u64, 254_u64),
        "ta:Z:chain"
    );

    format!("{}{}", first_half_gaf_line, second_half_gaf_line)
}

#[cfg(test)]
mod test {
    use std::path::PathBuf;

    use ab_poa::abpoa_wrapper::{AbpoaAligner, AbpoaAlignmentResult};
    use gfa::gfa::GFA;
    use gfa::parser::GFAParser;
    use handlegraph::handle::Edge;
    use handlegraph::hashgraph::HashGraph;
    use handlegraph::mutablehandlegraph::MutableHandleGraph;
    use handlegraph::pathgraph::PathHandleGraph;

    use crate::align::align_with_poa;
    use crate::align::find_nodes_edges_for_abpoa;
    use crate::align::generate_alignment;
    use crate::align::{find_range_chain, GAFAlignment};
    use crate::chain::{anchors_for_query, chain_anchors, Chain, Anchor, score_anchor};
    use crate::index::Index;
    use crate::io::QuerySequence;
    use crate::kmer::{SeqOrient, SeqPos};
    use ab_poa::abpoa::abpoa_align_sequence_to_subgraph;
    use handlegraph::handle::Handle;
    use crate::kmer::SeqOrient::Forward;
    //use rayon::prelude::*;

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
    fn anchors_found() {
        let graph = create_simple_graph();
        let index = Index::build(&graph, 3, 100, 100, None, None, false, None);

        let input_seq = QuerySequence::from_string(&String::from("ACTGCA"));
        let anchors = anchors_for_query(&index, &input_seq, true);
        assert_eq!(anchors.len(), 4); //NOT 2 because anchors can overlap

        let input_seq = QuerySequence::from_string(&String::from("AGAGC"));
        let anchors = anchors_for_query(&index, &input_seq, true);
        assert_eq!(anchors.len(), 3);
    }

    #[test]
    fn anchors_found_2() {
        let mut graph2 = HashGraph::new();
        let h1 = graph2.create_handle("AAAAAAAAAAA".as_bytes(), 1);
        let h2 = graph2.create_handle("C".as_bytes(), 2);
        let h3 = graph2.create_handle("G".as_bytes(), 3);
        let h4 = graph2.create_handle("TTTTTTTTTTTT".as_bytes(), 4);

        graph2.create_edge(&Edge(h1, h2));
        graph2.create_edge(&Edge(h1, h3));
        graph2.create_edge(&Edge(h2, h4));
        graph2.create_edge(&Edge(h3, h4));

        let kmer_length: u64 = 11;
        let index = Index::build(&graph2, kmer_length, 100, 100, None, None, false, None);

        let query = String::from("AAAAACTTTTTT");
        assert_eq!(query.len() as u64, kmer_length+1);
        let input_seq = QuerySequence::from_string(&query);
        let anchors = anchors_for_query(&index, &input_seq, true);
        assert_eq!(anchors.len(), 2);
        println!("Anchors are: {:#?}", anchors);
    }

    #[test]
    fn anchors_found_single_node() {
        let mut graph2 = HashGraph::new();
        let h1 = graph2.create_handle("AAATTAAA".as_bytes(), 1);

        let kmer_length: u64 = 3;
        let index = Index::build(&graph2, kmer_length, 100, 100, None, None, false, None);

        let query = String::from("AAAAAA");
        let input_seq = QuerySequence::from_string(&query);
        let mut anchors = anchors_for_query(&index, &input_seq, true);
        println!("Anchors are: {:#?}", anchors);
        let chains = chain_anchors(&mut anchors, kmer_length, 100, 100, 0,0.5, 10.0, 100.0, &input_seq);
        //println!("Chains are: {:#?}", chains);
    }


    #[test]
    fn test_simple_anchors() {
        let mut graph = HashGraph::new();
        graph.create_handle("ACT".as_bytes(), 1);
        let index = Index::build(&graph, 3, 100, 100, None, None, false, None);
        let query: String = String::from("ACT");
        let input_seq = QuerySequence::from_string(&query);

        let anchors = anchors_for_query(&index, &input_seq, false);
        assert_eq!(anchors.len(), 1);

        let anchor = anchors.get(0).unwrap();
        assert_eq!(anchor.query_begin, 0);
        assert_eq!(anchor.query_end, 3);
        assert_eq!(anchor.target_begin.position, 0);
        assert_eq!(anchor.target_begin.orient, SeqOrient::Forward);
        assert_eq!(anchor.target_end.position, 3);
        assert_eq!(anchor.target_end.orient, SeqOrient::Forward);
    }

    #[test]
    fn test_simple_anchors_reverse() {
        let mut graph = HashGraph::new();
        let h1 = graph.create_handle("AAA".as_bytes(), 1);
        let h2 = graph.create_handle("CCC".as_bytes(), 2);
        let h3 = graph.create_handle("GGG".as_bytes(), 3);
        let h4 = graph.create_handle("AAA".as_bytes(), 4);

        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h1, h3));
        graph.create_edge(&Edge(h2, h4));
        graph.create_edge(&Edge(h3, h4));

        let index = Index::build(&graph, 3, 100, 100, None, None, false, None);

        let input_seq = QuerySequence::from_string(&String::from("TTT"));
        let anchors = anchors_for_query(&index, &input_seq, false);
        assert_eq!(anchors.len(), 2);

        let handle_anchor_0 = index.handle_from_seqpos(&anchors.get(0).unwrap().target_begin);
        assert_eq!(u64::from(handle_anchor_0.id()), 4);
        assert!(handle_anchor_0.is_reverse());
        assert_eq!(
            handle_anchor_0,
            index.handle_from_seqpos(&anchors.get(0).unwrap().get_end_seqpos_inclusive())
        );

        let handle_anchor_1 = index.handle_from_seqpos(&anchors.get(1).unwrap().target_begin);
        assert_eq!(u64::from(handle_anchor_1.id()), 1);
        assert!(handle_anchor_1.is_reverse());
        assert_eq!(
            handle_anchor_1,
            index.handle_from_seqpos(&anchors.get(1).unwrap().get_end_seqpos_inclusive())
        );
    }

    #[test]
    fn test_simple_anchors_reverse_2() {
        let mut graph = HashGraph::new();
        let h1 = graph.create_handle("AAA".as_bytes(), 1);
        let h2 = graph.create_handle("CCC".as_bytes(), 2);
        let h3 = graph.create_handle("GGG".as_bytes(), 3);
        let h4 = graph.create_handle("AAA".as_bytes(), 4);

        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h1, h3));
        graph.create_edge(&Edge(h2, h4));
        graph.create_edge(&Edge(h3, h4));

        let index = Index::build(&graph, 9, 100, 100, None, None, false, None);

        let input_seq = QuerySequence::from_string(&String::from("TTTCCCTTT"));
        let anchors = anchors_for_query(&index, &input_seq, false);
        assert_eq!(anchors.len(), 1);

        let handle_anchor_start = index.handle_from_seqpos(&anchors.get(0).unwrap().target_begin);
        assert_eq!(u64::from(handle_anchor_start.id()), 4);
        assert!(handle_anchor_start.is_reverse());

        let handle_anchor_end =
            index.handle_from_seqpos(&anchors.get(0).unwrap().get_end_seqpos_inclusive());
        assert_eq!(u64::from(handle_anchor_end.id()), 1);
        assert!(handle_anchor_end.is_reverse());
    }

    #[test]
    fn test_anchors() {
        let graph = create_simple_graph();
        let index = Index::build(&graph, 3, 100, 100, None, None, false, None);
        let input_seq = QuerySequence::from_string(&String::from("ACTGCA"));
        let anchors = anchors_for_query(&index, &input_seq, false);

        // There are at least 4 3-mers in the query, there can be more because
        // a kmer can appear in multiple positions
        assert!(anchors.len() >= 4);
    }

    #[test]
    fn test_no_anchors() {
        let graph = create_simple_graph();
        let index = Index::build(&graph, 3, 100, 100, None, None, false, None);
        let input_seq = QuerySequence::from_string(&String::from("AAATTT"));
        let anchors = anchors_for_query(&index, &input_seq, false);
        assert_eq!(anchors.len(), 0)
    }

    #[test]
    fn test_no_anchors_2() {
        let graph = create_simple_graph();
        let index = Index::build(&graph, 3, 100, 100, None, None, false, None);
        let input_seq = QuerySequence::from_string(&String::from(""));
        let anchors = anchors_for_query(&index, &input_seq, false);
        assert_eq!(anchors.len(), 0)
    }

    #[test]
    fn test_chains() {
        let graph = create_simple_graph();
        let index = Index::build(&graph, 3, 100, 100, None, None, false, None);
        let input_seq = QuerySequence::from_string(&String::from("ACTGCA"));
        let mut anchors = anchors_for_query(&index, &input_seq, false);

        //println!("Anchors: {:#?}", anchors);
        println!("Len: {:#?}", anchors.len());
        let chains: Vec<Chain> = chain_anchors(
            &mut anchors,
            index.kmer_length,
            50,
            1000,
            1,
            0.5f64,
            0.1f64,
            60.0f64,
            &QuerySequence::new(),
        );
        //assert!(chains.is_empty());
        //println!("Chains: {:#?}", chains);
        //println!("Chains length: {}", chains.len());
    }

    #[test]
    fn test_chains_2() {
        // Create HashGraph from GFA
        let parser = GFAParser::new();
        let gfa: GFA<usize, ()> = parser
            .parse_file(&PathBuf::from("./test/test.gfa"))
            .unwrap();
        let graph = HashGraph::from_gfa(&gfa);

        let index = Index::build(&graph, 11, 100, 100, None, None, false, None);
        let input_seq = QuerySequence::from_string(&String::from(index.seq_fwd.clone()));
        let mut anchors = anchors_for_query(&index, &input_seq, false);

        //println!("Anchors len: {}", anchors.len());
        //println!("Anchors: {:#?}", anchors);
        assert!(!anchors.is_empty());

        let chains: Vec<Chain> = chain_anchors(
            &mut anchors,
            index.kmer_length,
            50,
            1000,
            2,
            0.5f64,
            0.1f64,
            60.0f64,
            &QuerySequence::new(),
        );
        assert!(!chains.is_empty());
        //println!("Chains_2: {:#?}", chains);
        //println!("Chains_2 length: {}", chains.len());
    }

    /*
    #[test]
    fn test_no_chains() {
        let graph = create_simple_graph();
        let index = Index::build(&graph, 3, 100, 100, None, None, false);
        let input_seq = QuerySequence::from_string(&String::from("AAATTT"));
        let mut anchors = anchors_for_query(&index, &input_seq);
        let chains: Vec<Chain> = chain_anchors(
            &mut anchors,
            index.kmer_length,
            50,
            1000,
            3,
            0.5f64,
            0.1f64,
            60.0f64,
            &QuerySequence::new(),
        );
        assert!(chains.is_empty());
    }
     */

    #[test]
    fn test_score_anchors() {
        let a = Anchor {
            id: 36,
            query_begin: 35,
            query_end: 46,
            target_begin: SeqPos {
                orient: Forward,
                position: 3907,
            },
            target_end: SeqPos {
                orient: Forward,
                position: 3918,
            },
            max_chain_score: 31.397,
            best_predecessor_id: None,
        };

        let b = Anchor {
            id: 51,
            query_begin: 49,
            query_end: 60,
            target_begin: SeqPos {
                orient: Forward,
                position: 3906,
            },
            target_end: SeqPos {
                orient: Forward,
                position: 3918,
            },
            max_chain_score: 49.0,
            best_predecessor_id: None,
        };

        assert_eq!(-f64::MAX, score_anchor(&a,&b, &11, &100))
    }
}
