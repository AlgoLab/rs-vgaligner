use std::collections::VecDeque;

use bstr::ByteVec;
use bv::*;
use gfa::gfa::Orientation;
use handlegraph::handle::{Direction, Handle, NodeId};
use handlegraph::handlegraph::HandleGraph;
use handlegraph::hashgraph::HashGraph;
use rayon::iter::ParallelIterator;
use rayon::prelude::{IntoParallelIterator, IntoParallelRefIterator, ParallelSliceMut};
use serde::{Deserialize, Serialize};

/// Additional data about each node, to be used together with seq_bv
#[derive(Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct NodeRef {
    /// Starting position on the forward
    pub(crate) seq_idx: u64,
    /// Starting position on the edge vector
    pub(crate) edge_idx: u64,
    /// Number of incoming edges to a node
    pub(crate) edges_to_node: u64,
}

/// Find the length of the sequence encoded by the graph
pub fn find_graph_seq_length(graph: &HashGraph) -> u64 {
    let graph_handles: Vec<Handle> = graph.handles_iter().collect();
    graph_handles
        .into_par_iter()
        .map(|h| graph.sequence(h).len() as u64)
        .sum()
}

/// Find the forward sequence encoded in a (not necessarily partially ordered) graph
/// This function uses a BFS visit, and retrieves the forward sequence from each
/// node it visits, then concatenates them into a string.
/// It also requires a new bivector to be passed as input, this will
/// contain the starting positions of each node in the forward (and reverse) sequence.
pub fn find_forward_sequence_bfs(graph: &HashGraph, seq_bv: &mut BitVec) -> String {
    let mut forward: String = String::new();
    let mut bv_pos: u64 = 0;

    // Create queue
    // NOTE: this is a Queue based BFS implementation, this was done
    // in order not to get a stack overflow
    let mut q: VecDeque<NodeId> = VecDeque::new();
    let mut visited_nodes_list: Vec<NodeId> = Vec::new();

    // Insert first value
    q.push_back(graph.min_id);

    while !q.is_empty() {
        let curr_node_id = q.pop_front().unwrap();
        let curr_node = graph.get_node(&curr_node_id).unwrap();
        forward.push_str(curr_node.sequence.to_string().as_str());

        // set bitvector
        seq_bv.set(bv_pos, true); //end marker
        bv_pos += curr_node.sequence.to_string().len() as u64;

        let current_handle = Handle::new(curr_node_id, Orientation::Forward);
        for neighbor in graph.handle_edges_iter(current_handle, Direction::Right) {
            if !visited_nodes_list.contains(&neighbor.id()) {
                // Add to visited nodes
                visited_nodes_list.push(neighbor.id());

                // Add new node to queue
                q.push_back(neighbor.id());
            }
        }
    }

    // Set end marker for bitvector
    seq_bv.set(bv_pos, true);

    forward
}

/// Find the forward sequence in a partially ordered graph, by following
/// the order of the handles. Also computes the bitvector representing
/// the node start positions in the forward, and the Node References.
pub fn find_forward_sequence(
    graph: &HashGraph,
    seq_bv: &mut BitVec,
    node_ref: &mut Vec<NodeRef>,
    graph_edges: &mut Vec<Handle>,
    n_edges: &mut u64,
) -> String {
    let mut forward: String = String::new();
    let mut bv_pos: u64 = 0;
    //let mut edge_pos: u64 = 0;

    // Sort both handles and edges since this is a PO
    let mut graph_handles: Vec<Handle> = graph.handles_iter().collect();
    graph_handles.par_sort();

    // Iterate over the handles
    for current_handle in graph_handles {
        // Find the sequence
        let curr_seq = graph.sequence(current_handle).into_string_lossy();
        forward.push_str(curr_seq.as_str());

        // Find nodes at the left of current node
        let left_edges_handles: Vec<Handle> = graph
            .handle_edges_iter(current_handle, Direction::Left)
            .collect();
        let left_edges_count = left_edges_handles.par_iter().count();

        // Create node_ref
        let curr_node_ref = NodeRef {
            seq_idx: bv_pos,
            edge_idx: *n_edges,
            edges_to_node: left_edges_count as u64,
        };
        node_ref.push(curr_node_ref);

        // Add left handles to graph_edges
        graph_edges.extend(left_edges_handles);
        *n_edges += left_edges_count as u64;

        // Find nodes at the right of current node
        let right_edges_handles: Vec<Handle> = graph
            .handle_edges_iter(current_handle, Direction::Right)
            .collect();
        let right_edges_count = right_edges_handles.par_iter().count();
        // Add right handles to graph_edges
        graph_edges.extend(right_edges_handles);
        *n_edges += right_edges_count as u64;

        // Set bitvector
        seq_bv.set(bv_pos, true); //end marker
        bv_pos += curr_seq.to_string().len() as u64;
    }

    // Set end marker for bitvector
    seq_bv.set(bv_pos, true);

    // Set end marker for edges
    let last_node_ref = NodeRef {
        seq_idx: bv_pos,
        edge_idx: *n_edges,
        edges_to_node: 0,
    };
    node_ref.push(last_node_ref);

    forward
}

/*

#[test]
fn test_simple_rank_1_a() {
    assert_eq!(get_bv_rank(&bit_vec![false]), vec![0])
}

#[test]
fn test_simple_rank_1_b() {
    assert_eq!(get_bv_rank(&bit_vec![true]), vec![1])
}

#[test]
fn test_simple_rank_2() {
    assert_eq!(get_bv_rank(&bit_vec![false, true, false]), vec![0, 1, 1])
}

#[test]
fn test_simple_rank_3() {
    assert_eq!(
        get_bv_rank(&bit_vec![true, false, true, false]),
        vec![1, 1, 2, 2]
    )
}
 */
