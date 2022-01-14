use core::cmp;
use std::ops::Range;

use ab_poa::abpoa_wrapper::{AbpoaAligner, AbpoaAlignmentResult};
use handlegraph::handle::Handle;
use itertools::Itertools;
use rayon::prelude::*;

use crate::chain::Chain;
use crate::index::Index;

use log::{info, warn};
use std::env;
use std::time::Instant;

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
        .map(|chain| match chain.is_placeholder {
            false => obtain_base_level_alignment(index, chain),
            true => GAFAlignment::from_placeholder_chain(chain),
        })
        .collect();
    alignments.par_sort_by(|a, b| b.path_length.cmp(&a.path_length));
    //println!("Alignments: {:#?}", alignments);
    alignments.first().cloned().unwrap()
}

/// Get the POA alignment starting from a [chain].
pub(crate) fn obtain_base_level_alignment(index: &Index, chain: &Chain) -> GAFAlignment {
    // Find the range of node ids involved in the alignment

    let start_range = Instant::now();
    let po_range = find_range_chain(index, chain);
    info!(
        "Finding the PO range took: {} ms",
        start_range.elapsed().as_millis()
    );
    //println!("Graph range: {:#?}", po_range);

    // Find nodes and edges
    let start_find_graph = Instant::now();
    let (nodes, edges) = find_nodes_edges_for_abpoa(&index, &po_range);
    info!(
        "Finding nodes and edges took: {} ms",
        start_find_graph.elapsed().as_millis()
    );
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
        //let start_alignment = Instant::now();
        result = align_with_poa(&nodes_str, &edges, chain.query.seq.as_str());
        /*
        info!(
            "Performing the alignment took: {} ms",
            start_alignment.elapsed().as_millis()
        );
         */
    }
    let start_GAF = Instant::now();
    let alignment: GAFAlignment = generate_alignment(
        index,
        chain,
        &result,
        &po_range,
        //&subquery_range,
        &(0 as u64..chain.query.seq.len() as u64),
        chain.query.seq.len(),
    );
    info!(
        "Generating the GAF took: {} ms",
        start_GAF.elapsed().as_millis()
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
    query_name: Option<String>,
    query_length: Option<u64>,
    query_start: Option<u64>,
    query_end: Option<u64>,
    strand: Option<char>,
    path_matching: Option<String>,
    path_length: Option<u64>,
    path_start: Option<u64>,
    path_end: Option<u64>,
    residue: Option<u64>,
    alignment_block_length: Option<u64>,
    mapping_quality: Option<u64>,
    notes: Option<String>,
}
impl GAFAlignment {
    pub(crate) fn from_chain(chain: &Chain, index: &Index) -> Self {
        assert_eq!(chain.is_placeholder, false);

        // Create a path matching (list of oriented nodes) for the chain.
        // Since we don't have an actual path on the graph, we just take
        // the nodes where the anchors are.
        let chain_path_matching: Vec<String> = chain
            .anchors
            .par_iter()
            .map(|anchor| {
                // Find first and last node, and output it in the path matching
                let first_handle = index.handle_from_seqpos(&anchor.target_begin);
                let last_handle = index.handle_from_seqpos(&anchor.target_end);

                // Add the orientation
                let first_node_str: String = match first_handle.is_reverse() {
                    false => ">".to_string() + u64::from(first_handle.id()).to_string().as_str(),
                    true => "<".to_string() + u64::from(first_handle.id()).to_string().as_str(),
                };
                let last_node_str: String = match last_handle.is_reverse() {
                    false => ">".to_string() + u64::from(last_handle.id()).to_string().as_str(),
                    true => "<".to_string() + u64::from(last_handle.id()).to_string().as_str(),
                };

                // If both ends of the anchor are on the same node, return it
                // only once
                match first_node_str.eq(&last_node_str) {
                    false => first_node_str + last_node_str.as_str(),
                    true => first_node_str,
                }
            })
            .collect();

        GAFAlignment {
            query_name: Some(chain.query.name.clone()),
            query_length: Some(chain.query.seq.len() as u64),
            query_start: Some(chain.anchors.front().unwrap().query_begin),
            query_end: Some(chain.anchors.back().unwrap().query_end),
            strand: Some('+'),
            path_matching: Some(chain_path_matching.join(",")),
            path_length: Some(0),
            path_start: Some(0),
            path_end: Some(0),
            residue: Some(0),
            alignment_block_length: Some(0),
            mapping_quality: Some(cmp::min(chain.mapping_quality as u64, 254_u64)),
            notes: Some("ta:Z:chain".to_string()),
        }
    }

    pub(crate) fn from_placeholder_chain(chain: &Chain) -> Self {
        assert_eq!(chain.is_placeholder, true);
        GAFAlignment {
            query_name: Some(chain.query.name.clone()),
            query_length: Some(chain.query.seq.len() as u64),
            query_start: None,
            query_end: None,
            strand: None,
            path_matching: None,
            path_length: None,
            path_start: None,
            path_end: None,
            residue: None,
            alignment_block_length: None,
            mapping_quality: Some(0),
            notes: None,
        }
    }

    pub fn to_string(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            match &self.query_name {
                Some(x) => x.as_str(),
                _ => &"*",
            },
            match &self.query_length {
                Some(x) => x.to_string(),
                _ => String::from("*"),
            },
            match &self.query_start {
                Some(x) => x.to_string(),
                _ => String::from("*"),
            },
            match &self.query_end {
                Some(x) => x.to_string(),
                _ => String::from("*"),
            },
            match &self.strand {
                Some(x) => x.to_string(),
                _ => String::from("*"),
            },
            match &self.path_matching {
                Some(x) => x.as_str(),
                _ => &"*",
            },
            match &self.path_length {
                Some(x) => x.to_string(),
                _ => String::from("*"),
            },
            match &self.path_start {
                Some(x) => x.to_string(),
                _ => String::from("*"),
            },
            match &self.path_end {
                Some(x) => x.to_string(),
                _ => String::from("*"),
            },
            match &self.residue {
                Some(x) => x.to_string(),
                _ => String::from("*"),
            },
            match &self.alignment_block_length {
                Some(x) => x.to_string(),
                _ => String::from("*"),
            },
            match &self.mapping_quality {
                Some(x) => x.to_string(),
                _ => String::from("*"),
            },
            match &self.notes {
                Some(x) => x.as_str(),
                _ => &"*",
            }
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

    let start_abPOA_build = Instant::now();
    aligner.add_nodes_edges(nodes, edges);
    info!(
        "Building the abPOA graph took: {} ms",
        start_abPOA_build.elapsed().as_millis()
    );

    let start_abPOA_alignment = Instant::now();
    let res = aligner.align_sequence(query);
    info!(
        "Alignment with abPOA took: {} ms",
        start_abPOA_alignment.elapsed().as_millis()
    );

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
        query_name: Some(chain.query.name.clone()),
        query_length: Some(subquery_range.end - subquery_range.start),
        query_start: Some(subquery_range.start),
        query_end: Some(subquery_range.end),
        strand: Some('+'),
        path_matching: Some(alignment_path_string.iter().join("")),
        path_length: Some(result.abpoa_nodes.len() as u64),
        path_start: Some(0),                                 //og_path_start,
        path_end: Some(result.abpoa_nodes.len() as u64 - 1), //og_path_end,
        residue: Some(0),
        alignment_block_length: Some(result.n_aligned_bases as u64),
        mapping_quality: Some(255), //result.best_score as u64,
        notes: Some("as:i:-30".to_string() + " " + result.cs_string.as_str()), //+ " " + "cg:Z:" + result.cigar.as_str()),
    }
}

#[cfg(test)]
mod test {
    use crate::align::GAFAlignment;
    use crate::chain::Chain;
    use crate::io::QuerySequence;

    #[test]
    fn test_to_string_placeholder() {
        // Simulate input read
        let read = QuerySequence::from_name_and_string("Read1", "AAACTA");

        // Create placeholder chain from that read
        let mut c = Chain::new_placeholder();
        c.query = read.clone();

        // Test that the alignment is printed out correctly
        let alignment = GAFAlignment::from_placeholder_chain(&c);
        let expected_string = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            read.name,
            read.seq.len(),
            "*",
            "*",
            "*",
            "*",
            "*",
            "*",
            "*",
            "*",
            "*",
            0,
            "*"
        );
        assert_eq!(alignment.to_string(), expected_string);
    }
}
