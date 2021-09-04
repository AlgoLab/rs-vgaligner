use crate::chain::Chain;
use crate::index::Index;
use crate::io::InputSequence;
use ab_poa::abpoa_wrapper::{AbpoaAligner, AbpoaAlignmentResult};
use bstr::ByteVec;
use handlegraph::handle::Handle;
use handlegraph::handlegraph::HandleGraph;
use handlegraph::hashgraph::HashGraph;
use itertools::Itertools;
use rayon::prelude::*;
use std::ops::Range;
use substring::Substring;
use core::cmp;

/// Find the leftmost starting node(handle) and rightmost node(handle) for each chain
pub fn extract_graph_po_range(index: &Index, chains: &Vec<Chain>) -> Vec<Range<u64>> {
    let start_end_vec: Vec<Range<u64>> = chains
        .par_iter()
        .map(|c| find_range_single_chain(index, c))
        .collect();
    start_end_vec
}
fn find_range_single_chain(index: &Index, chain: &Chain) -> Range<u64> {
    let min_node = chain
        .anchors
        .par_iter()
        .map(|a| index.get_node_from_pos(a.target_begin as usize, a.target_begin_orient))
        .min()
        .unwrap();
    let max_node = chain
        .anchors
        .par_iter()
        .map(|a| index.get_node_from_pos(a.target_end as usize, a.target_end_orient))
        .max()
        .unwrap();
    min_node..max_node
}
pub fn extract_query_subsequence(
    index: &Index,
    chains: &Vec<Chain>,
    _query: &InputSequence,
) -> Vec<String> {
    let subseqs: Vec<String> = chains
        .par_iter()
        .map(|c| {
            extract_single_subsequence(
                index,
                c.target_begin as usize,
                c.target_begin_orient,
                c.target_end as usize,
                c.target_end_orient,
            )
        })
        .collect();
    subseqs
}
fn extract_single_subsequence(
    index: &Index,
    begin: usize,
    begin_orient: bool,
    end: usize,
    end_orient: bool,
) -> String {
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
#[derive(Debug)]
pub(crate) struct GAFAlignment {
    query_name: String,
    query_length: u64,
    query_start: u64,
    query_end: u64,
    strand: char,
    path_matching: String,
    path_length: u64,
    path_start: u64,
    path_end: u64,
    residue: u64,
    alignment_block_length: u64,
    mapping_quality: u64,
    notes: String
}
impl GAFAlignment {
    pub(crate) fn from_chain(chain: &Chain, query: &InputSequence) -> Self {
        GAFAlignment {
            query_name: query.name.clone(),
            query_length: query.seq.len() as u64,
            query_start: chain.anchors.front().unwrap().query_begin,
            query_end: chain.anchors.back().unwrap().query_end,
            strand: '+',
            path_matching: " ".to_string(),
            path_length: 0,
            path_start: 0,
            path_end: 0,
            residue: 0,
            alignment_block_length:  0,
            mapping_quality: cmp::min(chain.mapping_quality as u64, 254_u64),
            notes: "ta:Z:chain".to_string()
        }
    }

    pub fn to_string(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            self.query_name,
            self.query_length,
            self.query_start,
            self.query_end,
            self.strand,
            self.path_matching,
            self.path_length,
            self.path_start,
            self.path_end,
            self.residue,
            self.alignment_block_length,
            self.mapping_quality,
            self.notes
        )
    }
}
pub(crate) fn obtain_base_level_alignment(
    index: &Index,
    chain: &Chain,
    query: &InputSequence,
) -> GAFAlignment {
    // Find the range implied by the chain (on the graph!!!)
    let po_range = find_range_single_chain(index, chain);

    // Find nodes and edges
    let (nodes, edges) = find_nodes_edges_2(&index, &po_range);

    // TODO: possibly avoid this
    let nodes_str: Vec<&str> = nodes.par_iter().map(|x| &x as &str).collect();
    let edges_usize: Vec<(usize, usize)> = edges
        .into_par_iter()
        .map(|x| (x.0 as usize, x.1 as usize))
        .collect();

    /* -- this is wrong, the whole query should be aligned --
    // Find subquery implied by the chain
    let subquery_range = chain.find_query_start_end();
    let subquery = query
        .to_string()
        .substring(subquery_range.start as usize, subquery_range.end as usize)
        .to_string();
     */

    // Align with abpoa
    let mut result = AbpoaAlignmentResult::new();
    unsafe {
        result = align_with_poa(&nodes_str, &edges_usize, query.seq.as_str());
    }
    let alignment: GAFAlignment = generate_alignment(chain, &result);
    alignment
}

/*
fn find_nodes_edges(graph: &HashGraph, po_range: &Range<u64>) -> (Vec<String>, Vec<(u64, u64)>) {
    let mut handles: Vec<Handle> = graph.handles_iter().collect();
    handles.par_sort();
    let mut node_handles: Vec<Handle> = Vec::new();
    for val in po_range.start..po_range.end {
        node_handles.push(*handles.get((val - 1) as usize).unwrap());
    }
    let mut node_seqs: Vec<String> = Vec::new();
    for handle in &node_handles {
        node_seqs.push(graph.sequence(*handle).into_string_lossy());
    }
    let mut edges: Vec<(u64, u64)> = node_handles
        .iter()
        .combinations(2)
        .filter_map(|h| {
            let left = **h.get(0).unwrap();
            let right = **h.get(1).unwrap();
            if graph.has_edge(left, right) {
                Some((u64::from(left), u64::from(left)))
            } else {
                None
            }
        })
        .collect();
    // Since in abpoa nodes start from 0, the edges have to be "normalized".
    // In practice this means subtracting from each each edge the start of the po_range,
    // e.g. po_range: [4,7), edges: [(4,5), (5,6), (6,7)] -> norm_edges [(0,1), (1,2), (2, 3)]
    let norm_edges: Vec<(u64, u64)> = edges
        .into_par_iter()
        .map(|x| {
            let start = x.0;
            let end = x.1;
            (start - po_range.start, end - po_range.start)
        })
        .collect();
    (node_seqs, norm_edges)
}
*/

/// Finds the nodes and edges involved in the alignment
fn find_nodes_edges_2(index: &Index, po_range: &Range<u64>) -> (Vec<String>, Vec<(u64, u64)>) {
    let seqs: Vec<String> = po_range.clone().into_par_iter().map(|x| {
        let handle_ref = index.node_ref.get((x-1) as usize).unwrap();
        let next_handle_ref = index.node_ref.get((x) as usize).unwrap();
        index.seq_fwd.substring(handle_ref.seq_idx as usize, next_handle_ref.seq_idx as usize).to_string()
    }).collect();

    // Generate edges
    // TODO: check indices
    let edges: Vec<(u64, u64)> = po_range.clone().into_par_iter().flat_map(|x|{
        let handle_ref = index.node_ref.get((x-1) as usize).unwrap();
        let next_handle_ref = index.node_ref.get((x) as usize).unwrap();
        (handle_ref.edge_idx..next_handle_ref.edge_idx).into_par_iter().map(move |edge| {
            // Since in abpoa nodes start from 0, the edges have to be "normalized".
            // In practice this means subtracting from each each edge the start of the po_range,
            // e.g. po_range: [4,7), edges: [(4,5), (5,6), (6,7)] -> norm_edges [(0,1), (1,2), (2, 3)]
            (x - po_range.start, edge - po_range.start)
        })
    }).collect();

    (seqs, edges)
}

/// Use rs-abpoa to align the sequences to the graph
unsafe fn align_with_poa(
    nodes: &Vec<&str>,
    edges: &Vec<(usize, usize)>,
    query: &str,
) -> AbpoaAlignmentResult {
    let mut aligner = AbpoaAligner::new_with_example_params();
    aligner.add_nodes_edges(nodes, edges);
    let res = aligner.align_sequence(query);
    res
}
// TODO: complete this
// The final alignment should be made up by:
// - perfect matches (= anchors) where applicable
// - poa result (cigar) otherwise
fn generate_alignment(chain: &Chain, result: &AbpoaAlignmentResult) -> GAFAlignment {
    GAFAlignment {
        query_name: "".to_string(),
        query_length: 0,
        query_start: 0,
        query_end: 0,
        strand: '+',
        path_matching: "".to_string(),
        path_length: 0,
        path_start: 0,
        path_end: 0,
        residue: 0,
        alignment_block_length: 0,
        mapping_quality: 0,
        notes: "".to_string()
    }
}
