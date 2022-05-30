use core::cmp;
use std::collections::HashMap;
use std::ops::Range;

use ab_poa::abpoa_wrapper::{AbpoaAligner, AbpoaAlignmentMode, AbpoaAlignmentResult};
use handlegraph::handle::{Direction, Edge, Handle};
use itertools::Itertools;
//use rayon::prelude::*;

use crate::chain::Chain;
use crate::index::Index;

use crate::kmer::SeqOrient;
use crate::validate::{create_subgraph_GFA, export_GFA};
use ab_poa::abpoa::rand;
use handlegraph::hashgraph::{HashGraph, PathId};
use handlegraph::pathgraph::PathHandleGraph;
use log::{info, warn};
use std::env;
use std::time::Instant;
use handlegraph::handlegraph::HandleGraph;
use handlegraph::mutablehandlegraph::MutableHandleGraph;
use rspoa::api::{align_global_no_gap, align_local_no_gap};
use rspoa::gaf_output::GAFStruct;
use seal::pair::NeedlemanWunsch;

/// Get all the alignments from the [query_chains], but only return the best one.
pub fn best_alignment_for_query(
    index: &Index,
    query_chains: &Vec<Chain>,
    align_best_n: u64,
    graph: &HashGraph,
    export_subgraphs: bool,
) -> GAFAlignment {
    //println!("Query: {}, Chains: {:#?}", query_chains.get(0).unwrap().query.seq, query_chains);
    let mut alignments: Vec<GAFAlignment> = query_chains
        .iter()
        .take(cmp::min(align_best_n as usize, query_chains.len()))
        .map(|chain| match chain.is_placeholder {
            false => obtain_base_level_alignment(index, chain, graph, export_subgraphs),
            true => GAFAlignment::from_placeholder_chain(chain),
        })
        .collect();

    alignments.sort_by(|a, b| b.path_length.cmp(&a.path_length));
    //println!("Alignments: {:#?}", alignments);
    alignments.first().cloned().unwrap()
}

/// Get the POA alignment starting from a [chain].
pub(crate) fn obtain_base_level_alignment(
    index: &Index,
    chain: &Chain,
    graph: &HashGraph,
    export_subgraph: bool,
) -> GAFAlignment {
    info!("Start alignment of {}!", chain.query.name);

    // Find the range of node ids involved in the alignment
    let start_range = Instant::now();
    let po_range = find_range_chain(index, chain);
    info!(
        "Finding the PO range took: {} ms",
        start_range.elapsed().as_millis()
    );
    println!(
        "Graph range (before): {:?}",
        po_range
            .handles
            .iter()
            .map(|x| x.unpack_number())
            .collect::<Vec<u64>>()
    );
    let extended_range = extend_range_chain_2(index, chain, po_range);
    println!(
        "Graph range (after): {:?}",
        extended_range
            .handles
            .iter()
            .map(|x| x.unpack_number())
            .collect::<Vec<u64>>()
    );

    // Find nodes and edges
    let start_find_graph = Instant::now();
    let (nodes, edges) = find_nodes_edges_for_abpoa(&index, &extended_range);
    info!(
        "Finding nodes and edges took: {} ms",
        start_find_graph.elapsed().as_millis()
    );
    // TODO: possibly avoid this
    let nodes_str: Vec<&str> = nodes.iter().map(|x| &x as &str).collect();

    // The obtained subgraph can be exported as a GAF, if the user
    // chooses to do so (mostly for debug purposes)
    if export_subgraph {
        let paths = get_subgraph_paths(graph, &extended_range);
        let subgraph_as_gfa = create_subgraph_GFA(&nodes_str, &edges, &paths);
        export_GFA(
            subgraph_as_gfa,
            format!("{}-subgraph-{}.gfa", chain.query.name, chain.anchors.len()),
        )
        .unwrap();

        info!(
            "GFA exported in {}",
            format!(
                "subgraphs/{}-subgraph-{}.gfa",
                chain.query.name, chain.score
            )
        );
    }

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

    // Align with rspoa
    let subgraph = generate_subgraph_hashgraph(nodes_str, edges);
    for h in subgraph.handles_iter().sorted() {
        println!("Node is: {}", h.unpack_number());
        println!("Neighbors are: {:?}", subgraph.handle_edges_iter(h, Direction::Right).map(|x| x.unpack_number()).collect::<Vec<u64>>());
    }

    /*
    let res_gaf = align_global_no_gap(
        &chain.query.seq.to_string(),
        &subgraph,
        None, None, Some(1.0f32));
     */
    let res_gaf = align_local_no_gap(
        &chain.query.seq.to_string(),
        &subgraph,
        None,
        None);
    //println!("Res GAF: {}", res_gaf.to_string());

    let alignment = GAFAlignment::from_rspoa_alignment(res_gaf, chain, extended_range);

    /*
    let alignment = GAFAlignment {
        query_name: Some(chain.query.name.clone()),
        query_length: None,
        query_start: None,
        query_end: None,
        strand: None,
        path_matching: None,
        path_length: None,
        path_start: None,
        path_end: None,
        residue: None,
        alignment_block_length: None,
        mapping_quality: None,
        notes: None
    };
     */

    /*
    // Align with abpoa
    let result: AbpoaAlignmentResult;
    unsafe {
        //result = align_with_poa(&nodes_str, &edges, subquery.as_str());
        //let start_alignment = Instant::now();
        //result = align_with_poa(&nodes_str, &edges, chain.query.seq.as_str());
        /*
        info!(
            "Performing the alignment took: {} ms",
            start_alignment.elapsed().as_millis()
        );
         */

        let alignment_mode = match nodes_str.len() {
            1 => AbpoaAlignmentMode::Global,
            _ => AbpoaAlignmentMode::Local
        };

        let mode_as_str = match alignment_mode {
            AbpoaAlignmentMode::Global => "Global",
            AbpoaAlignmentMode::Local => "Local"
        };

        info!(
        "For read {} calling abPOA with: \n nodes: vec!{:?} \n edges: vec!{:?} \n query: \"{}\" \n mode: {}",
        chain.query.name, nodes_str, edges, chain.query.seq.as_str(), mode_as_str
        );

        result = AbpoaAligner::create_align_safe(&nodes_str, &edges, chain.query.seq.as_str(), alignment_mode);
    }

    //Test seal
    //let strategy = NeedlemanWunsch::new(1, -1, -1, -1);

    println!("Abpoa result abpoa-nodes {:?}", result.abpoa_nodes);
    println!("Abpoa result graph-nodes {:?}", result.graph_nodes);

    let start_GAF = Instant::now();
    let alignment: GAFAlignment = generate_alignment(
        index,
        chain,
        &result,
        &extended_range,
        //&subquery_range,
        &(0 as u64..chain.query.seq.len() as u64),
        chain.query.seq.len(),
    );
    info!(
        "Generating the GAF took: {} ms",
        start_GAF.elapsed().as_millis()
    );
     */

    alignment
}

pub fn generate_subgraph_hashgraph(nodes_str: Vec<&str>, edges: Vec<(usize, usize)>) -> HashGraph {
    let mut subgraph = HashGraph::new();

    let handles: Vec<Handle> = nodes_str.iter().map(|seq| subgraph.append_handle(seq.as_bytes())).collect();

    for edge in edges {
        let left_handle = handles.get(edge.0).unwrap();
        let right_handle = handles.get(edge.1).unwrap();
        subgraph.create_edge(&Edge(*left_handle, *right_handle))
    }

    subgraph
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
        .iter()
        .map(|a| index.handle_from_seqpos(&a.target_begin))
        .collect();
    let end_handles: Vec<Handle> = chain
        .anchors
        .iter()
        // Here a method is called (instead of simply getting the attribute)
        // because the end position is stored as non-inclusive, but I want the
        // last included position
        .map(|a| index.handle_from_seqpos(&a.get_end_seqpos_inclusive()))
        .collect();

    //println!("Start ids: {:?}", start_handles.iter().map(|x| x.unpack_number()).collect::<Vec<u64>>());
    //println!("End ids: {:?}", end_handles.iter().map(|x| x.unpack_number()).collect::<Vec<u64>>());

    let all_handles: Vec<Handle> = start_handles
        .into_iter()
        .chain(end_handles.into_iter())
        .collect();
    let min_handle: Handle = *all_handles.iter().min().unwrap();
    let max_handle: Handle = *all_handles.iter().max().unwrap();

    //println!("ids: {:?}", all_handles.iter().map(|x| x.unpack_number()).collect::<Vec<u64>>());
    //println!("Min id: {}", min_handle.unpack_number());
    //println!("Max id: {}", max_handle.unpack_number());

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
        po_range_handles.sort();

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

pub fn extend_range_chain(
    index: &Index,
    chain: &Chain,
    old_range: OrientedGraphRange,
) -> OrientedGraphRange {
    let mut extended_handles: Vec<Handle> = old_range.handles.clone();

    println!(
        "Chain start (query): {}, chain ends (query): {}",
        chain.anchors.front().unwrap().query_begin,
        chain.anchors.back().unwrap().query_end
    );

    let mut prefix_diff: u64 = chain.anchors.front().unwrap().query_begin;
    println!("Prefix diff is: {}", prefix_diff);

    // Only if necessary, extend to the left
    if prefix_diff > 0 {
        let first_handle = old_range.get_first_handle();

        let mut left_handles: Vec<(usize, Handle)> = index
            .incoming_edges_from_handle(first_handle)
            .into_iter()
            .map(|handle| (prefix_diff as usize, handle))
            .collect();

        while !left_handles.is_empty() {
            let mut next_handles: Vec<(usize, Handle)> = vec![];

            for (prefix_left, curr_handle) in left_handles {
                println!(
                    "Extending node {} with prefix_left: {}",
                    curr_handle.unpack_number(),
                    prefix_left
                );
                // First, push the handle (since the prefix_diff is > 0)
                extended_handles.push(curr_handle);

                // Then, determine if this handle has to be extended further
                let curr_handle_seq = index.seq_from_handle(&curr_handle);
                if curr_handle_seq.len() < prefix_left {
                    let remaining_length = prefix_left - curr_handle_seq.len();

                    // Find the left neighbours of the curr_handle
                    let new_left_handles: Vec<(usize, Handle)> = index
                        .incoming_edges_from_handle(&curr_handle)
                        .into_iter()
                        .map(|handle| (remaining_length as usize, handle))
                        .collect();

                    next_handles.extend(new_left_handles.into_iter());
                }

                // Otherwise, the curr_handle does not get re-added, and is therefore ignored
                // in the next iterations
            }

            left_handles = next_handles;
        }
    }

    let mut suffix_diff: u64 =
        chain.query.seq.len() as u64 - chain.anchors.back().unwrap().query_end;
    println!("Suffix diff is: {}", suffix_diff);

    // Then, if necessary, extend to the right
    if suffix_diff > 0 {
        let last_handle = old_range.get_last_handle();

        let mut right_handles: Vec<(usize, Handle)> = index
            .outgoing_edges_from_handle(last_handle)
            .into_iter()
            .map(|handle| (suffix_diff as usize, handle))
            .collect();

        while !right_handles.is_empty() {
            let mut next_handles: Vec<(usize, Handle)> = vec![];

            for (suffix_left, curr_handle) in right_handles {
                println!(
                    "Extending node {} with suffix_left: {}",
                    curr_handle.unpack_number(),
                    suffix_left
                );
                // First, push the handle (since the prefix_diff is > 0)
                extended_handles.push(curr_handle);

                // Then, determine if this handle has to be extended further
                let curr_handle_seq = index.seq_from_handle(&curr_handle);
                if curr_handle_seq.len() < suffix_left {
                    let remaining_length = suffix_left - curr_handle_seq.len();

                    // Find the left neighbours of the curr_handle
                    let new_right_handles: Vec<(usize, Handle)> = index
                        .outgoing_edges_from_handle(&curr_handle)
                        .into_iter()
                        .map(|handle| (remaining_length as usize, handle))
                        .collect();

                    next_handles.extend(new_right_handles.into_iter());
                }

                // Otherwise, the curr_handle does not get re-added, and is therefore ignored
                // in the next iterations
            }

            right_handles = next_handles;
        }
    }

    extended_handles.sort();
    extended_handles.dedup();

    OrientedGraphRange {
        orient: old_range.orient,
        handles: extended_handles,
    }
}

pub fn extend_range_chain_2(
    index: &Index,
    chain: &Chain,
    old_range: OrientedGraphRange,
) -> OrientedGraphRange {
    let mut extended_handles: Vec<Handle> = old_range.handles.clone();

    println!(
        "Chain start (query): {}, chain ends (query): {}",
        chain.anchors.front().unwrap().query_begin,
        chain.anchors.back().unwrap().query_end
    );

    let mut prefix_diff: u64 = chain.anchors.front().unwrap().query_begin;
    println!("Prefix diff is: {}", prefix_diff);

    let first_handle = old_range.get_first_handle();
    let first_anchor = chain.anchors.front().unwrap();
    let start_prefix_on_node =
        first_anchor.target_begin.position - index.get_bv_select(first_handle.unpack_number());
    if start_prefix_on_node < prefix_diff {
        prefix_diff -= start_prefix_on_node;
    } else {
        prefix_diff = 0;
    }
    println!("New prefix diff is: {}", prefix_diff);

    // Only if necessary, extend to the left
    if prefix_diff > 0 {
        let mut left_handles: Vec<(usize, Handle)> = index
            .incoming_edges_from_handle(first_handle)
            .into_iter()
            .map(|handle| (prefix_diff as usize, handle))
            .collect();

        while !left_handles.is_empty() {
            let mut next_handles: Vec<(usize, Handle)> = vec![];

            for (prefix_left, curr_handle) in left_handles {
                println!(
                    "Extending node {} with prefix_left: {}",
                    curr_handle.unpack_number(),
                    prefix_left
                );
                // First, push the handle (since the prefix_diff is > 0)
                extended_handles.push(curr_handle);

                // Then, determine if this handle has to be extended further
                let curr_handle_seq = index.seq_from_handle(&curr_handle);
                if curr_handle_seq.len() < prefix_left {
                    let remaining_length = prefix_left - curr_handle_seq.len();

                    // Find the left neighbours of the curr_handle
                    let new_left_handles: Vec<(usize, Handle)> = index
                        .incoming_edges_from_handle(&curr_handle)
                        .into_iter()
                        .map(|handle| (remaining_length as usize, handle))
                        .collect();

                    next_handles.extend(new_left_handles.into_iter());
                }

                // Otherwise, the curr_handle does not get re-added, and is therefore ignored
                // in the next iterations
            }

            left_handles = next_handles;
        }
    }

    let mut suffix_diff: u64 =
        chain.query.seq.len() as u64 - chain.anchors.back().unwrap().query_end;
    println!("Suffix diff is: {}", suffix_diff);

    let last_handle = old_range.get_last_handle();
    let last_anchor = chain.anchors.back().unwrap();
    let end_suffix_on_node = index.get_bv_select(last_handle.unpack_number() + 1)
        - 1
        - last_anchor.get_end_seqpos_inclusive().position;
    println!(
        "End suffix on node: {} (next_node_start-1: {}, anchor_end: {})",
        end_suffix_on_node,
        index.get_bv_select(first_handle.unpack_number() + 1) - 1,
        last_anchor.get_end_seqpos_inclusive().position
    );
    if end_suffix_on_node > suffix_diff {
        suffix_diff = 0;
    } else {
        suffix_diff -= end_suffix_on_node;
    }
    println!("New suffix diff is: {}", suffix_diff);

    // Then, if necessary, extend to the right
    if suffix_diff > 0 {
        let mut right_handles: Vec<(usize, Handle)> = index
            .outgoing_edges_from_handle(last_handle)
            .into_iter()
            .map(|handle| (suffix_diff as usize, handle))
            .collect();

        while !right_handles.is_empty() {
            let mut next_handles: Vec<(usize, Handle)> = vec![];

            for (suffix_left, curr_handle) in right_handles {
                println!(
                    "Extending node {} with suffix_left: {}",
                    curr_handle.unpack_number(),
                    suffix_left
                );
                // First, push the handle (since the prefix_diff is > 0)
                extended_handles.push(curr_handle);

                // Then, determine if this handle has to be extended further
                let curr_handle_seq = index.seq_from_handle(&curr_handle);
                if curr_handle_seq.len() < suffix_left {
                    let remaining_length = suffix_left - curr_handle_seq.len();

                    // Find the left neighbours of the curr_handle
                    let new_right_handles: Vec<(usize, Handle)> = index
                        .outgoing_edges_from_handle(&curr_handle)
                        .into_iter()
                        .map(|handle| (remaining_length as usize, handle))
                        .collect();

                    next_handles.extend(new_right_handles.into_iter());
                }

                // Otherwise, the curr_handle does not get re-added, and is therefore ignored
                // in the next iterations
            }

            right_handles = next_handles;
        }
    }

    extended_handles.sort();
    extended_handles.dedup();

    OrientedGraphRange {
        orient: old_range.orient,
        handles: extended_handles,
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
        .iter()
        .map(|h| index.seq_from_handle(h))
        .collect();

    // Get the edges between the handles in the po_range
    // In abpoa the nodes are given as 0-based tuples (starting_node, ending_node).
    // For example in order to create this graph: "ACG" -> "AAA"
    // the following instruction will be used:
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
                    let starting_pos = range_handles.iter().position(|x| *x == handle).unwrap();
                    let ending_pos = range_handles
                        .iter()
                        .position(|x| *x == target_handle)
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
    pub query_name: Option<String>,
    pub query_length: Option<u64>,
    pub query_start: Option<u64>,
    pub query_end: Option<u64>,
    pub strand: Option<char>,
    pub path_matching: Option<String>,
    pub path_length: Option<u64>,
    pub path_start: Option<u64>,
    pub path_end: Option<u64>,
    pub residue: Option<u64>,
    pub alignment_block_length: Option<u64>,
    pub mapping_quality: Option<u64>,
    pub notes: Option<String>,
}
impl GAFAlignment {
    pub(crate) fn from_chain(chain: &Chain, index: &Index) -> Self {
        assert_eq!(chain.is_placeholder, false);

        // Create a path matching (list of oriented nodes) for the chain.
        // Since we don't have an actual path on the graph, we just take
        // the nodes where the anchors are.
        /*
        let chain_path_matching: Vec<String> = chain
            .anchors
            .iter()
            .map(|anchor| {
                // Find first and last node, and output it in the path matching
                let first_handle = index.handle_from_seqpos(&anchor.target_begin);
                let last_handle = index.handle_from_seqpos(&anchor.get_end_seqpos_inclusive());

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
            */

        let first_chain_node: Option<u64> = None;
        let last_chain_node: Option<u64> = None;
        //info!("Anchors: {:#?}", chain.anchors);
        let mut chain_path_matching: Vec<String> = chain
            .anchors
            .iter()
            .map(|anchor| {
                /*
                 // Find first and last node, and output it in the path matching
                 let first_handle = index.handle_from_seqpos(&anchor.target_begin);
                 let last_handle = index.handle_from_seqpos(&anchor.get_end_seqpos_inclusive());

                 let first_node_start = index.get_bv_select(u64::from(first_handle.id()));
                 let first_node_offset = anchor.target_begin.position - first_node_start as u64;

                 let last_node_start = index.get_bv_select(u64::from(last_handle.id()));
                 let last_node_offset = anchor.target_end.position - last_node_start as u64;

                /*
                 let first_handle_rank = first_handle.unpack_number() - 1;
                 let first_node_start = index
                     .node_ref
                     .get(first_handle_rank as usize)
                     .unwrap()
                     .seq_idx;
                 let first_node_offset = anchor.target_begin.position - first_node_start;

                 let last_handle_rank = last_handle.unpack_number() - 1;
                 let last_node_start = index
                     .node_ref
                     .get(last_handle_rank as usize)
                     .unwrap()
                     .seq_idx;
                 let last_node_offset = anchor.target_end.position - last_node_start;
                  */

                 // Add the orientation
                 let first_node_str: String = match first_handle.is_reverse() {
                     false => {
                         ">".to_string()
                             + u64::from(first_handle.id()).to_string().as_str()
                             + ":"
                             + first_node_offset.to_string().as_str()
                     }
                     true => "<".to_string() + u64::from(first_handle.id()).to_string().as_str(),
                 };
                 let last_node_str: String = match last_handle.is_reverse() {
                     false => {
                         ">".to_string()
                             + u64::from(last_handle.id()).to_string().as_str()
                             + ":"
                             + last_node_offset.to_string().as_str()
                     }
                     true => "<".to_string() + u64::from(last_handle.id()).to_string().as_str(),
                 };

                 format!("({},{}),", first_node_str, last_node_str)
                  */

                let graph_pos = anchor.get_detailed_graph_position(index);

                let first_node_str: String = match graph_pos.start_orient {
                    SeqOrient::Forward => {
                        ">".to_string()
                            + u64::from(graph_pos.start_node).to_string().as_str()
                            + ":"
                            + graph_pos.start_offset_from_node.to_string().as_str()
                    }
                    SeqOrient::Reverse => {
                        "<".to_string()
                            + u64::from(graph_pos.start_node).to_string().as_str()
                            + ":"
                            + graph_pos.start_offset_from_node.to_string().as_str()
                    }
                };

                let last_node_str: String = match graph_pos.end_orient {
                    SeqOrient::Forward => {
                        ">".to_string()
                            + u64::from(graph_pos.end_node).to_string().as_str()
                            + ":"
                            + graph_pos.end_offset_from_node.to_string().as_str()
                    }
                    SeqOrient::Reverse => {
                        "<".to_string()
                            + u64::from(graph_pos.end_node).to_string().as_str()
                            + ":"
                            + graph_pos.end_offset_from_node.to_string().as_str()
                    }
                };

                format!("({},{}),", first_node_str, last_node_str)
            })
            .collect();

        GAFAlignment {
            query_name: Some(chain.query.name.clone()),
            query_length: Some(chain.query.seq.len() as u64),
            query_start: Some(chain.anchors.front().unwrap().query_begin),
            query_end: Some(chain.anchors.back().unwrap().query_end),
            strand: Some('+'),
            path_matching: Some(chain_path_matching.join("")),
            path_length: Some(0),
            path_start: Some(0),
            path_end: Some(0),
            residue: Some(0),
            alignment_block_length: Some(0),
            mapping_quality: Some(cmp::min(chain.mapping_quality as u64, 254_u64)),
            notes: Some(format!(
                "{},{}",
                "ta:Z:chain".to_string(),
                format!("n_anchors: {}", chain.anchors.len())
            )),
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

    pub(crate) fn from_rspoa_alignment(rspoa_alignment: GAFStruct, chain: &Chain, graph_range: OrientedGraphRange) -> Self {

        println!("Range handles are: {:?}", graph_range.handles);
        println!("Path matching nodes are: {:?}", rspoa_alignment.path);

        let og_graph_handles: Vec<Handle> =
            rspoa_alignment
                .path
                .iter()
                .map(|rspoa_node| {
                    graph_range.handles.get(*rspoa_node-1).unwrap().clone()
                })
                .collect();

        let alignment_path_string: Vec<String> = og_graph_handles
            .iter()
            .map(|handle| match handle.is_reverse() {
                false => ">".to_string() + u64::from(handle.id()).to_string().as_str(),
                true => "<".to_string() + u64::from(handle.id()).to_string().as_str(),
            })
            .collect();

        GAFAlignment {
            query_name: Some(chain.query.name.clone()),
            query_length: Some(chain.query.seq.len() as u64),
            query_start: Some(rspoa_alignment.query_start as u64),
            query_end: Some(rspoa_alignment.query_end as u64),
            strand: Some(rspoa_alignment.strand),
            path_matching: Some(alignment_path_string.join("")),
            path_length: Some(rspoa_alignment.path_length as u64),
            path_start: Some(rspoa_alignment.path_start as u64),
            path_end: Some(rspoa_alignment.path_end as u64),
            residue: Some(rspoa_alignment.residue_matches_number as u64),
            alignment_block_length: Some(0),
            mapping_quality: Some(255),
            notes: Some(rspoa_alignment.comments),
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
        .iter()
        .map(|x| *po_range.handles.get(*x).unwrap())
        .collect();
    //println!("Og graph aligment path: {:#?}", og_graph_alignment_path);
    //println!("Og graph nodes: {:#?}", og_graph_alignment_path.par_iter().map(|x| u64::from(x.id())).collect::<Vec<u64>>());

    let alignment_path_string: Vec<String> = og_graph_alignment_path
        .iter()
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
        //path_start: Some(0),                                 //og_path_start,
        //path_end: Some(result.abpoa_nodes.len() as u64 - 1), //og_path_end,
        path_start: Some(result.aln_start_offset as u64),
        path_end: Some(result.aln_end_offset as u64),
        residue: Some(0),
        alignment_block_length: Some(result.n_aligned_bases as u64),
        mapping_quality: Some(255), //result.best_score as u64,
        notes: Some(
            "as:i:-30".to_string()
                + " "
                + result.cs_string.as_str()
                + ",cg:Z:"
                + result.cigar.as_str(),
        ),
    }
}

pub fn get_subgraph_paths(
    graph: &HashGraph,
    po_range: &OrientedGraphRange,
) -> HashMap<PathId, Vec<u64>> {
    let mut subgraph_paths: HashMap<PathId, Vec<u64>> = HashMap::new();

    let min_in_range = po_range.handles.iter().min().unwrap().unpack_number();

    for path_id in graph.paths_iter() {
        let mut nodes_in_path: Vec<u64> = vec![];
        let curr_path = graph.get_path(path_id).unwrap();
        for handle in &curr_path.nodes {
            if po_range.handles.contains(&handle) {
                nodes_in_path.push(handle.unpack_number() - min_in_range + 1);
            }
        }
        subgraph_paths.insert(*path_id, nodes_in_path);
    }
    subgraph_paths
}

#[cfg(test)]
mod test {
    use crate::align::{get_subgraph_paths, GAFAlignment, OrientedGraphRange, RangeOrient};
    use crate::chain::Chain;
    use crate::io::QuerySequence;
    use gfa::gfa::{Orientation, GFA};
    use gfa::parser::GFAParser;
    use handlegraph::handle::Handle;
    use handlegraph::hashgraph::HashGraph;
    use std::hash::Hash;
    use std::path::PathBuf;

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

    #[test]
    fn get_graph_paths() {
        // Create HashGraph from GFA
        let parser = GFAParser::new();
        let gfa: GFA<usize, ()> = parser
            .parse_file(&PathBuf::from("./test/test.gfa"))
            .unwrap();
        let graph = HashGraph::from_gfa(&gfa);

        let oriented_graph_range = OrientedGraphRange {
            orient: RangeOrient::Forward,
            handles: vec![
                Handle::new(graph.min_id, Orientation::Forward),
                Handle::new(graph.max_id, Orientation::Forward),
            ],
        };

        println!(
            "Get subgraph path: {:#?}",
            get_subgraph_paths(&graph, &oriented_graph_range)
        );
    }
}
