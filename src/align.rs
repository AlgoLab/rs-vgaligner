use core::cmp;
use std::ops::Range;

use ab_poa::abpoa_wrapper::{AbpoaAligner, AbpoaAlignmentResult};
use handlegraph::handle::Handle;
use itertools::Itertools;
use rayon::prelude::*;

use crate::chain::Chain;
use crate::index::Index;

/// Get all the alignments from the [query_chains], but only return the best one.
pub fn best_alignment_for_query(
    index: &Index,
    query_chains: &Vec<Chain>,
    align_best_n: u64,
) -> GAFAlignment {
    //println!("Query: {}, Chains: {:#?}", query_chains.get(0).unwrap().query.seq, query_chains);
    let mut alignments: Vec<GAFAlignment> = query_chains
        .par_iter()
        .take(cmp::min(align_best_n as usize, query_chains.len()))
        .map(|chain| obtain_base_level_alignment(index, chain))
        .collect();
    alignments.par_sort_by(|a, b| b.path_length.cmp(&a.path_length));
    //println!("Alignments: {:#?}", alignments);
    alignments.first().cloned().unwrap()
}

/// Get the POA alignment starting from a [chain].
pub(crate) fn obtain_base_level_alignment(index: &Index, chain: &Chain) -> GAFAlignment {
    // Find the range of node ids involved in the alignment
    let po_range = find_range_chain(index, chain);
    //println!("Graph range: {:#?}", po_range);

    // Find nodes and edges
    let (nodes, edges) = find_nodes_edges_for_abpoa(&index, &po_range);
    // TODO: possibly avoid this
    let nodes_str: Vec<&str> = nodes.par_iter().map(|x| &x as &str).collect();
    //println!("Seqs: {:#?}\n, edges: {:#?}\n, Query: {:#?}", nodes_str, edges, chain.query.seq.to_string());

    //println!("Query seq: {:#?}", chain.query.seq);
    //println!("Seqs: {:?}, len {}, nested-len {}", nodes_str, nodes_str.len(), nodes_str.par_iter().map(|x| x.len()).sum::<usize>());
    //println!("Edges: {:?}, len {}", edges, edges.len());
    //println!("Chain: {:#?}", chain);
    //println!("po_range: {:#?}", po_range);

    // Find subquery implied by the chain
    /*
    let subquery_range = chain.find_query_start_end();
    let subquery = chain
        .query
        .seq
        .to_string()
        .substring(subquery_range.start as usize, subquery_range.end as usize)
        .to_string();
     */
    //println!("Subquery is: {:#?}", subquery);

    // Align with abpoa
    let result: AbpoaAlignmentResult;
    unsafe {
        //result = align_with_poa(&nodes_str, &edges, subquery.as_str());
        result = align_with_poa(&nodes_str, &edges, chain.query.seq.as_str());
    }
    let alignment: GAFAlignment = generate_alignment(
        index,
        chain,
        &result,
        &po_range,
        //&subquery_range,
        &(0 as u64..chain.query.seq.len() as u64),
        chain.query.seq.len(),
    );

    alignment
}

/// Represents the orientation of a range of nodes (= a subgraph but without edges) in a graph.
#[derive(Debug)]
pub enum RangeOrient {
    Forward,
    Reverse,
    Both,
}
/// Represents an oriented range of nodes of a graph.
#[derive(Debug)]
pub struct OrientedGraphRange {
    pub orient: RangeOrient,
    pub handles: Vec<Handle>,
}
impl OrientedGraphRange {
    pub fn get_first_handle(&self) -> &Handle {
        self.handles.get(0).unwrap()
    }
    pub fn get_last_handle(&self) -> &Handle {
        self.handles.last().unwrap()
    }
}

/// Find the range of nodes implied by the chain.
pub fn find_range_chain(index: &Index, chain: &Chain) -> OrientedGraphRange {
    /* -- old version: this did not work correctly for reverse handles
    let mut min_handle: Handle = chain
        .anchors
        .par_iter()
        .map(|a| index.handle_from_seqpos(&a.target_begin))
        .min()
        .unwrap();
    let mut max_handle: Handle = chain
        .anchors
        .par_iter()
        // Here a method is called (instead of simply getting the attribute)
        // because the end position is stored as non-inclusive, but I want the
        // last included position
        .map(|a| index.handle_from_seqpos(&a.get_end_seqpos_inclusive()))
        .max()
        .unwrap();
     */

    let start_handles: Vec<Handle> = chain
        .anchors
        .par_iter()
        .map(|a| index.handle_from_seqpos(&a.target_begin))
        .collect();
    let end_handles: Vec<Handle> = chain
        .anchors
        .par_iter()
        // Here a method is called (instead of simply getting the attribute)
        // because the end position is stored as non-inclusive, but I want the
        // last included position
        .map(|a| index.handle_from_seqpos(&a.get_end_seqpos_inclusive()))
        .collect();
    let all_handles: Vec<Handle> = start_handles
        .into_par_iter()
        .chain(end_handles.into_par_iter())
        .collect();
    let min_handle: Handle = *all_handles.par_iter().min().unwrap();
    let max_handle: Handle = *all_handles.par_iter().max().unwrap();

    /*
    println!(
        "Min id: {}, handle: {:#?}",
        min_handle.id(),
        min_handle.as_integer()
    );
    println!(
        "Max id: {}, handle: {:#?}",
        max_handle.id(),
        max_handle.as_integer()
    );
     */

    /*
    // If on reverse strand the range min and max could be reversed
    if min_handle > max_handle {
        // But this should never happen in other cases...
        assert!(min_handle.is_reverse());
        assert!(max_handle.is_reverse());
        let temp = min_handle;
        min_handle = max_handle;
        max_handle = temp;
    }
     */

    let mut po_range_handles: Vec<Handle>;
    let orient: RangeOrient;

    /*
    println!(
        "Min handle: {:#?}, Max handle: {:#?}",
        min_handle, max_handle
    );
     */
    //println!("Range as handle is: {:#?}", min_handle..max_handle);
    /*
    println!(
        "Range as nodes is: {:#?}",
        u64::from(min_handle)..u64::from(max_handle)
    );
    println!("Noderef {:#?}", index.noderef_from_handle(&min_handle));
    println!("Handle: {:#?}, u64: {:#?}, Handle from u64: {:#?}", min_handle, u64::from(min_handle), Handle::from_integer(u64::from(min_handle)));
    */

    if !min_handle.is_reverse() && !max_handle.is_reverse() {
        po_range_handles = (u64::from(min_handle)..=u64::from(max_handle))
            .into_iter()
            .map(|x| Handle::from_integer(x * 2))
            .filter(|x| !x.is_reverse())
            .collect();
        orient = RangeOrient::Forward;
    } else if min_handle.is_reverse() && max_handle.is_reverse() {
        po_range_handles = (u64::from(min_handle)..=u64::from(max_handle))
            .into_iter()
            .map(|x| Handle::from_integer(x * 2 + 1))
            .filter(|x| x.is_reverse())
            .collect();
        orient = RangeOrient::Reverse;
    } else {
        let po_range_handles_fwd: Vec<Handle> = (u64::from(min_handle)..=u64::from(max_handle))
            .into_iter()
            .map(|x| Handle::from_integer(x * 2))
            .filter(|x| !x.is_reverse())
            .collect();
        let po_range_handles_rev: Vec<Handle> = (u64::from(min_handle)..=u64::from(max_handle))
            .into_iter()
            .map(|x| Handle::from_integer(x * 2 + 1))
            .filter(|x| x.is_reverse())
            .collect();
        po_range_handles = [po_range_handles_fwd, po_range_handles_rev].concat();
        po_range_handles.par_sort();

        orient = RangeOrient::Both
    }

    // If min_handle == max_handle (i.e. the chain is on a single Handle), due to how
    // ranges work in Rust, it would return an empty range. This if tries to fix that.
    // I previously tried to set the range as inclusive to the right (i.e. ..= ) but
    // it would lead to index.noderef to go out of bounds. (This is probably due to
    // upper bounds being non-inclusive throughout the program).
    if po_range_handles.is_empty() && min_handle == max_handle {
        po_range_handles.push(min_handle);
    }

    OrientedGraphRange {
        orient,
        handles: po_range_handles,
    }
}

/// Finds the nodes and edges involved in the alignment. These are returned in a format
/// compatible with rs-abpoa (https://github.com/HopedWall/rs-abPOA), see the
/// "Interacting with abPOA" section for further details.
pub(crate) fn find_nodes_edges_for_abpoa(
    index: &Index,
    po_range: &OrientedGraphRange,
) -> (Vec<String>, Vec<(usize, usize)>) {
    let range_handles = po_range.handles.clone();
    //println!("Range handles is: {:#?}", range_handles);

    // Get the seqs for the handles in the po_range
    // This is necessary because in abpoa the nodes are specified by a Vec of sequences
    let seqs: Vec<String> = range_handles
        .par_iter()
        .map(|h| index.seq_from_handle(h))
        .collect();

    // Get the edges between the handles in the po_range
    // In abpoa the nodes are given as 0-based tuples (starting_node, ending_node).
    // For example in order to create this graph: "ACG" -> "AAA"
    // the following instuction will be used:
    // aligner.add_nodes_edges(&vec!["ACG", "AAA"], &vec![(0, 1)]);
    let edges: Vec<(usize, usize)> = range_handles
        .clone()
        .into_iter()
        // For each handle in the po range
        .flat_map(|handle| {
            index
                // Find all the outgoing edges (this actually returns the handles the edges point to)
                .outgoing_edges_from_handle(&handle)
                .into_iter()
                // Only keep the ones in the po range
                .filter(|target_handle| range_handles.contains(target_handle))
                // And find their 0-based position in the po range
                .map(|target_handle| {
                    let starting_pos = range_handles
                        .par_iter()
                        .position_any(|x| *x == handle)
                        .unwrap();
                    let ending_pos = range_handles
                        .par_iter()
                        .position_any(|x| *x == target_handle)
                        .unwrap();
                    (starting_pos, ending_pos)
                }) //(handle.as_integer() - po_range.par_iter().min().unwrap().as_integer(), target_handle.as_integer() - po_range.par_iter().min().unwrap().as_integer()))
                .collect::<Vec<(usize, usize)>>()
        })
        .collect();

    // Remove eventual loops
    // Since nodes are ordered, it should be enough to keep only
    // the edges where the start nodes has a lower id than
    // the end node
    let edges_without_loops: Vec<(usize, usize)> = match po_range.orient {
        RangeOrient::Forward => edges.into_iter().filter(|edge| edge.0 < edge.1).collect(),
        RangeOrient::Reverse => edges.into_iter().filter(|edge| edge.1 < edge.0).collect(),
        RangeOrient::Both => edges,
    };

    (seqs, edges_without_loops)
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

    NOTE: There can also be additional notes at the end.
*/
/// Represents an alignment in GAF format.
#[derive(Debug, Clone)]
pub struct GAFAlignment {
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
    notes: String,
}
impl GAFAlignment {
    pub(crate) fn from_chain(chain: &Chain, index: &Index) -> Self {
        // Find first and last node, and output it in the path matching
        let first_chain_handle =
            index.handle_from_seqpos(&chain.anchors.front().unwrap().target_begin);
        let last_chain_handle = index.handle_from_seqpos(&chain.anchors.back().unwrap().target_end);

        let first_node_str: String = match first_chain_handle.is_reverse() {
            false => ">".to_string() + u64::from(first_chain_handle.id()).to_string().as_str(),
            true => "<".to_string() + u64::from(first_chain_handle.id()).to_string().as_str(),
        };
        let last_node_str: String = match last_chain_handle.is_reverse() {
            false => ">".to_string() + u64::from(last_chain_handle.id()).to_string().as_str(),
            true => "<".to_string() + u64::from(last_chain_handle.id()).to_string().as_str(),
        };

        GAFAlignment {
            query_name: chain.query.name.clone(),
            query_length: chain.query.seq.len() as u64,
            query_start: chain.anchors.front().unwrap().query_begin,
            query_end: chain.anchors.back().unwrap().query_end,
            strand: '+',
            path_matching: first_node_str + last_node_str.as_str(),
            path_length: 0,
            path_start: 0,
            path_end: 0,
            residue: 0,
            alignment_block_length: 0,
            mapping_quality: cmp::min(chain.mapping_quality as u64, 254_u64),
            notes: "ta:Z:chain".to_string(),
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
        .par_iter()
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

/// Use rs-abpoa to align the sequences to the graph.
pub(crate) unsafe fn align_with_poa(
    nodes: &Vec<&str>,
    edges: &Vec<(usize, usize)>,
    query: &str,
) -> AbpoaAlignmentResult {
    let mut aligner = AbpoaAligner::new_with_example_params();
    aligner.add_nodes_edges(nodes, edges);
    let res = aligner.align_sequence(query);
    res
}

/// Generate a GAF alignment from the result obtained from abPOA.
pub(crate) fn generate_alignment(
    index: &Index,
    chain: &Chain,
    result: &AbpoaAlignmentResult,
    po_range: &OrientedGraphRange,
    subquery_range: &Range<u64>,
    og_query_length: usize,
) -> GAFAlignment {
    //println!("CIGAR: {:#?}", result.cigar);

    // Get the nodes involved in the alignment from abPOA
    let mut alignment_path = result.graph_nodes.clone();
    //println!("Alignment path is: {:#?}\n", alignment_path);

    // Since abpoa deals with 1-base nodes, the abstraction simply converts
    // the abpoa_node_id to their position in the abstraction
    // i.e suppose the abpoa_nodes are 2,3,4,5; the abstraction nodes are [2] and [3,4,5];
    // result.graph_nodes will contain [0, 1, 1, 1], I want: [0,1]
    alignment_path.dedup();

    // Since alignment_path contains ids that are coherent with the nodes/edges
    // used during its creation, the nodes will be something like [0, 1, ... po_range.max].
    // However, this graph is a subgraph of our original graph, so the ids must be shifted
    // accordingly.
    let og_graph_alignment_path: Vec<Handle> = alignment_path
        .par_iter()
        .map(|x| *po_range.handles.get(*x).unwrap())
        .collect();
    //println!("Og graph aligment path: {:#?}", og_graph_alignment_path);
    //println!("Og graph nodes: {:#?}", og_graph_alignment_path.par_iter().map(|x| u64::from(x.id())).collect::<Vec<u64>>());

    let alignment_path_string: Vec<String> = og_graph_alignment_path
        .par_iter()
        .map(|handle| match handle.is_reverse() {
            false => ">".to_string() + u64::from(handle.id()).to_string().as_str(),
            true => "<".to_string() + u64::from(handle.id()).to_string().as_str(),
        })
        .collect();

    // Find path start and end
    /*
    let og_path_start_handle = og_graph_alignment_path.first().unwrap();
    let og_path_end_handle = og_graph_alignment_path.last().unwrap();
    let og_path_start = 0 +
            result.graph_nodes.iter().filter(|x| (**x as u64) == og_path_start_handle.as_integer()).count() as u64;
    let og_path_end =
            result.graph_nodes.iter().filter(|x| (**x as u64) == og_path_end_handle.as_integer()).count() as u64;
    */

    GAFAlignment {
        query_name: chain.query.name.clone(),
        query_length: subquery_range.end - subquery_range.start,
        query_start: subquery_range.start,
        query_end: subquery_range.end,
        strand: '+',
        path_matching: alignment_path_string.iter().join(""),
        path_length: result.abpoa_nodes.len() as u64,
        path_start: 0,                                 //og_path_start,
        path_end: result.abpoa_nodes.len() as u64 - 1, //og_path_end,
        residue: 0,
        alignment_block_length: result.n_aligned_bases as u64,
        mapping_quality: 255, //result.best_score as u64,
        notes: ("as:i:-30".to_string() + " " + result.cs_string.as_str()), //+ " " + "cg:Z:" + result.cigar.as_str()),
    }
}

#[cfg(test)]
mod test {

    #[test]
    fn test_find_node_edges() {}
}
