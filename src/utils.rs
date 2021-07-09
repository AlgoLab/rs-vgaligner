use std::collections::VecDeque;

use bv::*;
use gfa::gfa::Orientation;
use handlegraph::handle::{Direction, Edge, Handle, NodeId};
use handlegraph::handlegraph::HandleGraph;
use handlegraph::hashgraph::HashGraph;
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
    let mut total_length = 0;
    //for value in graph.handles_iter() {
    for value in graph.handles_iter() {
        let node = graph.get_node(&value.id()).unwrap();
        total_length += node.sequence.len() as u64;
    }
    total_length
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
    node_ref: &mut Vec<NodeRef>
) -> String {
    let mut forward: String = String::new();
    let mut bv_pos: u64 = 0;
    let mut edge_pos: u64 = 0;

    // Sort both handles and edges since this is a PO
    let mut graph_handles: Vec<Handle> = graph.handles_iter().collect();
    graph_handles.sort();

    let mut graph_edges: Vec<Edge> = graph.edges_iter().collect();
    graph_edges.sort();

    // Iterate over the handles
    for current_handle in graph_handles {
        // Find the sequence
        let curr_node_id = current_handle.id();
        let curr_node = graph.get_node(&curr_node_id).unwrap();
        forward.push_str(curr_node.sequence.to_string().as_str());

        // Create nodeRef
        let incoming_edges: u64 = graph
            .handle_edges_iter(current_handle, Direction::Left)
            .count() as u64;
        let curr_node_ref = NodeRef {
            seq_idx: bv_pos,
            edge_idx: edge_pos,
            edges_to_node: incoming_edges,
        };
        node_ref.push(curr_node_ref);

        // Advance the edge position so that it points to edges
        // relative to the next handle
        while edge_pos < graph_edges.len() as u64 {
            let curr_edge = graph_edges.get(edge_pos as usize).unwrap();
            if curr_edge.0 != current_handle {
                break;
            }
            edge_pos += 1;
        }

        // Set bitvector
        seq_bv.set(bv_pos, true); //end marker
        bv_pos += curr_node.sequence.to_string().len() as u64;
    }

    // Set end marker for bitvector
    seq_bv.set(bv_pos, true);

    forward
}

/// Return the rank vector of a given bitvector bv
/// The rank for an element bv[i] is defined as the number of 1s in [bv[0]...bv[i]]
/// We compute the rank for every element in bv, and put these ranks in a Vec
pub fn get_bv_rank(bv: &BitVec) -> Vec<u32> {
    let mut bv_rank = Vec::with_capacity(bv.len() as usize);

    let mut curr_rank = 0;
    for i in 0..bv.len() {
        if bv.get(i) == true {
            curr_rank += 1;
        }
        bv_rank.push(curr_rank);
    }

    bv_rank
}

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
