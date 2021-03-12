use bv::*;
use handlegraph::hashgraph::HashGraph;
use handlegraph::handle::{Handle, Direction, NodeId, Edge};
use gfa::gfa::Orientation;
use handlegraph::handlegraph::HandleGraph;
use std::collections::VecDeque;
use handlegraph::mutablehandlegraph::MutableHandleGraph;
use bstr::ByteSlice;

/// Finds the forward sequence encoded in a (not necessarily partially ordered) graph
/// This function uses a BFS visit, and retrieves the forward sequence from each
/// node it visits, then concatenates them into a string.
/// It also requires a new bivector to be passed as input, this will
/// contain the starting positions of each node in the forward (and reverse) sequence.
pub fn find_sequence(graph: &HashGraph, seq_bv: &mut BitVec) -> String {
    let mut forward : String = String::new();
    let mut bv_pos : u64 = 0;

    // Create queue
    // NOTE: this is a Queue based BFS implementation, this was done
    // in order not to get a stack overflow
    let mut q: VecDeque<NodeId> = VecDeque::new();
    let mut visited_nodes_list : Vec<NodeId> = Vec::new();

    // Insert first value
    q.push_back(graph.min_id);

    while !q.is_empty() {

        let curr_node_id = q.pop_front().unwrap();
        let curr_node = graph.get_node(&curr_node_id).unwrap();
        forward.push_str(curr_node.sequence.to_string().as_str());

        // set bitvector
        seq_bv.set(bv_pos, true); //end marker
        bv_pos += (curr_node.sequence.to_string().len() as u64);

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

/// Returns the rank vector of a given bitvector bv
/// The rank for an element bv[i] is defined as the number of 1s in [bv[0]...bv[i]]
/// We compute the rank for every element in bv, and put these ranks in a Vec
pub fn get_bv_rank(bv : &BitVec) -> Vec<u32> {
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
    assert_eq!(get_bv_rank(&bit_vec![false, true, false]), vec![0,1,1])
}

#[test]
fn test_simple_rank_3() {
    assert_eq!(get_bv_rank(&bit_vec![true, false, true, false]), vec![1,1,2,2])
}


