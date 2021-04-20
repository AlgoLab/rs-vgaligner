use crate::utils::NodeRef;
use ahash::AHashMap;
use ahash::CallHasher;
use ahash::RandomState;
use boomphf::*;
use bstr::ByteVec;
use bv::BitVec;
use handlegraph::handle::{Direction, Handle};
use handlegraph::handlegraph::HandleGraph;
use handlegraph::hashgraph::HashGraph;
use handlegraph::pathgraph::PathHandleGraph;
use std::cmp::min;
use std::collections::VecDeque;
use substring::Substring;

/// Struct that represents a kmer in the graph
#[derive(Debug, Clone)]
pub struct Kmer {
    /// The sequence of the kmer
    pub(crate) seq: String,
    /// The start position relative to the handle
    begin: u64,
    /// The end position relative to the handle
    end: u64,
    /// The handles where the kmer is
    pub(crate) handles: Vec<Handle>,
    /// The orientation of the handles
    handle_orient: bool,
    /// The number of forks the kmer has been involved in
    forks: u64,
}

impl PartialEq for Kmer {
    fn eq(&self, other: &Self) -> bool {
        self.seq == other.seq
    }
}

impl Kmer {
    /// This function allows for a kmer to be extended, used when the
    /// kmer is on more than one handle
    fn extend_kmer(&mut self, new_seq: String, new_handle: Handle) {
        self.seq.push_str(&new_seq);
        //self.end += new_seq.len() as u64;
        self.end = new_seq.len() as u64;
        self.handles.push(new_handle);
    }

    fn add_handle_to_complete(&mut self, new_handle: Handle) {
        self.handles.push(new_handle);
    }
}

/*pub fn generate_kmers_2(graph : &HashGraph, k : u64, edge_max : Option<u64>, degree_max : Option<u64>) -> Vec<Kmer> {
    let mut kmers : Vec<Kmer> = Vec::new();
    let mut curr_kmers_to_complete : Vec<Kmer> = Vec::new();

    let mut q: VecDeque<NodeId> = VecDeque::new();
    let mut visited_nodes : Vec<NodeId> = Vec::new();

    q.push_back(graph.min_id);
    while let Some(curr_node_id) = q.pop_front() {

        let current_handle = Handle::pack(curr_node_id, false);

        // First try extending the open kmers with the current node


        // Get current handle
        let node = graph.get_node(&curr_node_id).unwrap();
        let handle_seq = node.sequence.to_string();
        let handle_length = handle_seq.len() as u64;

        // Then try generating the kmers from the given node
        for i in 0..handle_length {
            let begin = i;
            let end = min(i+k, handle_length);
            let mut kmer = Kmer {
                seq : handle_seq.substring(begin as usize, end as usize).to_string(),
                begin,
                end,
                handle
            };

            if (kmer.seq.len() as u64) == k {
                if !kmers.contains(&kmer) {
                    kmers.push(kmer);
                }
            } else {
                curr_kmers_to_complete.push(kmer);
            }
        }


        for neighbor in graph.handle_edges_iter(current_handle, Direction::Right) {
            if !visited_nodes.contains(*neighbor.id()) {
                // Add neighbor id to queue
                q.push_back(neighbor.id());
            }
        }

        // Add to visited nodes
        visited_nodes.push(curr_node);
    }


    kmers
}*/

pub fn generate_kmers(
    graph: &HashGraph,
    k: u64,
    edge_max: Option<u64>,
    degree_max: Option<u64>,
) -> Vec<Kmer> {
    let mut complete_kmers: Vec<Kmer> = Vec::new();

    let mut sorted_graph_handles: Vec<Handle> = graph.handles_iter().collect();
    sorted_graph_handles.sort();

    for forward_handle in sorted_graph_handles {
        let mut handle: Handle;
        let mut orient: bool;

        for handle_orient in &[true, false] {
            match handle_orient {
                true => {
                    handle = forward_handle;
                    orient = true;
                }
                false => {
                    handle = forward_handle.flip();
                    orient = false;
                }
            }

            if let Some(degree_max) = degree_max {
                let mut curr_count: u64 = 0;
                graph
                    .handle_edges_iter(handle, Direction::Right)
                    .for_each(|_| {
                        curr_count += 1;
                    });
                if curr_count > degree_max {
                    // Skip current orientation
                    continue;
                }
            }

            // Get current node/handle
            let mut handle_seq = graph.sequence(handle).into_string_lossy();
            let mut handle_length = handle_seq.len() as u64;

            // Then try generating the kmers from the given node
            let mut incomplete_kmers: Vec<Kmer> = Vec::new();

            for i in 0..handle_length {
                let begin = i;
                let end = min(i + k, handle_length);
                let mut kmer = Kmer {
                    seq: handle_seq
                        .substring(begin as usize, end as usize)
                        .to_string(),
                    begin,
                    end,
                    handles: vec![handle],
                    handle_orient: orient,
                    forks: 0,
                };

                if (kmer.seq.len() as u64) == k {
                    if !complete_kmers.contains(&kmer) {
                        complete_kmers.push(kmer);
                    }
                } else {
                    let mut next_count: u64 = 0;
                    if edge_max.is_some() || degree_max.is_some() {
                        graph
                            .handle_edges_iter(handle, Direction::Right)
                            .for_each(|_| {
                                next_count += 1;
                            });
                    }

                    if (degree_max.is_none() && edge_max.is_none())
                        || (degree_max.is_some() && next_count < degree_max.unwrap())
                        || (edge_max.is_some() && kmer.forks < edge_max.unwrap())
                    {
                        for neighbor in graph.handle_edges_iter(handle, Direction::Right) {
                            let mut inc_kmer = kmer.clone();
                            inc_kmer.add_handle_to_complete(neighbor);

                            if next_count > 1 {
                                inc_kmer.forks += 1;
                            }

                            incomplete_kmers.push(inc_kmer);
                        }
                    }
                }
            }

            // Then complete all incomplete kmers
            while let Some(mut incomplete_kmer) = incomplete_kmers.pop() {
                handle = incomplete_kmer.handles.pop().unwrap();
                handle_seq = graph.sequence(handle).into_string_lossy();
                handle_length = handle_seq.len() as u64;

                let end = min(k - (incomplete_kmer.seq.len() as u64), handle_length);
                let str_to_add = handle_seq.substring(0, end as usize).to_string();
                incomplete_kmer.extend_kmer(str_to_add, handle);

                if (incomplete_kmer.seq.len() as u64) == k {
                    if !complete_kmers.contains(&incomplete_kmer) {
                        complete_kmers.push(incomplete_kmer);
                    }
                } else {
                    // NOTE: if there is no neighbor, the kmer does not get re-added
                    // to the incomplete ones, so that the external loop can end
                    for neighbor in graph.handle_edges_iter(handle, Direction::Right) {
                        let mut next_count: u64 = 0;
                        if edge_max.is_some() || degree_max.is_some() {
                            graph
                                .handle_edges_iter(handle, Direction::Right)
                                .for_each(|_| {
                                    next_count += 1;
                                });
                        }

                        if (degree_max.is_none() && edge_max.is_none())
                            || (degree_max.is_some() && next_count < degree_max.unwrap())
                            || (edge_max.is_some() && incomplete_kmer.forks < edge_max.unwrap())
                        {
                            let mut inc_kmer = incomplete_kmer.clone();
                            inc_kmer.add_handle_to_complete(neighbor);

                            if next_count > 1 {
                                inc_kmer.forks += 1;
                            }

                            incomplete_kmers.push(inc_kmer);
                        }
                    }
                }
            }

            assert!(incomplete_kmers.is_empty());
        }
    }

    complete_kmers
}

pub fn generate_kmers_linearly_forward(
    graph: &HashGraph,
    k: u64,
    edge_max: Option<u64>,
    degree_max: Option<u64>,
) -> Vec<Kmer> {
    let mut kmers: Vec<Kmer> = Vec::new();
    let mut prev_kmers_to_complete: VecDeque<Kmer> = VecDeque::new();

    //println!("Paths are: {:#?}", graph.paths.values());
    for path_id in graph.paths_iter() {
        let path = graph.get_path(path_id).unwrap();

        //println!("Path: {:#?}", path);
        prev_kmers_to_complete.clear();

        for handle in &path.nodes {
            let mut curr_kmers_to_complete: Vec<Kmer> = Vec::new();

            // Get current handle
            let handle_seq = graph.sequence(*handle).into_string_lossy();
            let handle_length = handle_seq.len() as u64;

            // First try completing previously uncompleted kmers
            // EXAMPLE: suppose k=4
            // node i-1 = AC, node i = CGT
            // incomplete kmer is AC (len 2), I want ACCG (len 4 = k)
            // also C + CGT
            while let Some(mut incomplete_kmer) = prev_kmers_to_complete.pop_front() {
                let end = min(k - (incomplete_kmer.seq.len() as u64), handle_length);
                let str_to_add = handle_seq.substring(0, end as usize).to_string();
                incomplete_kmer.extend_kmer(str_to_add, *handle);

                if (incomplete_kmer.seq.len() as u64) == k {
                    if !kmers.contains(&incomplete_kmer) {
                        kmers.push(incomplete_kmer);
                    }
                } else {
                    curr_kmers_to_complete.push(incomplete_kmer);
                }
            }

            // Then try to generate current node kmers
            for i in 0..handle_length {
                let begin = i;
                let end = min(i + k, handle_length);
                let mut kmer = Kmer {
                    seq: handle_seq
                        .substring(begin as usize, end as usize)
                        .to_string(),
                    begin,
                    end,
                    handles: vec![*handle],
                    handle_orient: true,
                    forks: 0,
                };

                if (kmer.seq.len() as u64) == k {
                    if !kmers.contains(&kmer) {
                        kmers.push(kmer);
                    }
                } else {
                    curr_kmers_to_complete.push(kmer);
                }
            }

            // Add kmers not completed in this iteration
            curr_kmers_to_complete
                .iter()
                .for_each(|inc_kmer| prev_kmers_to_complete.push_back(inc_kmer.clone()));
        }
    }

    kmers
}

pub fn generate_kmers_linearly_reverse(
    graph: &HashGraph,
    k: u64,
    edge_max: Option<u64>,
    degree_max: Option<u64>,
) -> Vec<Kmer> {
    let mut kmers: Vec<Kmer> = Vec::new();
    let mut prev_kmers_to_complete: VecDeque<Kmer> = VecDeque::new();

    //println!("Paths are: {:#?}", graph.paths.values());
    for path_id in graph.paths_iter() {
        let path = graph.get_path(path_id).unwrap();

        //println!("Path: {:#?}", path);
        prev_kmers_to_complete.clear();

        let mut reverse_order_handles = path.nodes.clone();
        reverse_order_handles.reverse();

        for handle_forward in reverse_order_handles {
            let handle = handle_forward.flip();

            let mut curr_kmers_to_complete: Vec<Kmer> = Vec::new();

            // Get current handle
            let handle_seq = graph.sequence(handle).into_string_lossy();
            let handle_length = handle_seq.len() as u64;

            // First try completing previously uncompleted kmers
            // EXAMPLE: suppose k=4
            // node i-1 = AC, node i = CGT
            // incomplete kmer is AC (len 2), I want ACCG (len 4 = k)
            // also C + CGT
            while let Some(mut incomplete_kmer) = prev_kmers_to_complete.pop_front() {
                let end = min(k - (incomplete_kmer.seq.len() as u64), handle_length);
                let str_to_add = handle_seq.substring(0, end as usize).to_string();
                incomplete_kmer.extend_kmer(str_to_add, handle);

                if (incomplete_kmer.seq.len() as u64) == k {
                    if !kmers.contains(&incomplete_kmer) {
                        kmers.push(incomplete_kmer);
                    }
                } else {
                    curr_kmers_to_complete.push(incomplete_kmer);
                }
            }

            // Then try to generate current node kmers
            for i in 0..handle_length {
                let begin = i;
                let end = min(i + k, handle_length);
                let mut kmer = Kmer {
                    seq: handle_seq
                        .substring(begin as usize, end as usize)
                        .to_string(),
                    begin,
                    end,
                    handles: vec![handle],
                    handle_orient: false,
                    forks: 0,
                };

                if (kmer.seq.len() as u64) == k {
                    if !kmers.contains(&kmer) {
                        kmers.push(kmer);
                    }
                } else {
                    curr_kmers_to_complete.push(kmer);
                }
            }

            // Add kmers not completed in this iteration
            curr_kmers_to_complete
                .iter()
                .for_each(|inc_kmer| prev_kmers_to_complete.push_back(inc_kmer.clone()));
        }
    }

    kmers
}

pub fn merge_kmers(kmers_fwd: Vec<Kmer>, kmers_rev: Vec<Kmer>) -> Vec<Kmer> {
    let mut kmers: Vec<Kmer> = kmers_fwd.clone();

    for kmer in kmers_rev {
        if !kmers.contains(&kmer) {
            kmers.push(kmer);
        }
    }

    kmers
}

pub fn generate_kmers_linearly(
    graph: &HashGraph,
    k: u64,
    edge_max: Option<u64>,
    degree_max: Option<u64>,
) -> Vec<Kmer> {
    let mut kmers: Vec<Kmer> = Vec::new();

    let forward_kmers = generate_kmers_linearly_forward(graph, k, edge_max, degree_max);
    let reverse_kmers = generate_kmers_linearly_reverse(graph, k, edge_max, degree_max);

    kmers = merge_kmers(forward_kmers, reverse_kmers);

    kmers
}

pub fn generate_kmers_linearly_2(
    graph: &HashGraph,
    k: u64,
    edge_max: Option<u64>,
    degree_max: Option<u64>,
) -> Vec<Kmer> {
    let mut complete_kmers: Vec<Kmer> = Vec::new();

    for path_id in graph.paths_iter() {
        let path = graph.get_path(path_id).unwrap();

        for forward_handle in &path.nodes {
            let mut handle: Handle;
            let mut orient: bool;

            for handle_orient in &[true, false] {
                match handle_orient {
                    true => {
                        handle = *forward_handle;
                        orient = true;
                    }
                    false => {
                        handle = forward_handle.flip();
                        orient = false;
                    }
                }

                if let Some(degree_max) = degree_max {
                    let mut curr_count: u64 = 0;
                    graph
                        .handle_edges_iter(handle, Direction::Right)
                        .for_each(|_| {
                            curr_count += 1;
                        });
                    if curr_count > degree_max {
                        // Skip current orientation
                        continue;
                    }
                }

                // Get current node/handle
                let mut handle_seq = graph.sequence(handle).into_string_lossy();
                let mut handle_length = handle_seq.len() as u64;

                // Then try generating the kmers from the given node
                let mut incomplete_kmers: Vec<Kmer> = Vec::new();

                for i in 0..handle_length {
                    let begin = i;
                    let end = min(i + k, handle_length);
                    let mut kmer = Kmer {
                        seq: handle_seq
                            .substring(begin as usize, end as usize)
                            .to_string(),
                        begin,
                        end,
                        handles: vec![handle],
                        handle_orient: orient,
                        forks: 0,
                    };

                    if (kmer.seq.len() as u64) == k {
                        if !complete_kmers.contains(&kmer) {
                            complete_kmers.push(kmer);
                        }
                    } else {
                        let mut next_count: u64 = 0;
                        if edge_max.is_some() || degree_max.is_some() {
                            graph
                                .handle_edges_iter(handle, Direction::Right)
                                .for_each(|_| {
                                    next_count += 1;
                                });
                        }

                        if (degree_max.is_none() && edge_max.is_none())
                            || (degree_max.is_some() && next_count < degree_max.unwrap())
                            || (edge_max.is_some() && kmer.forks < edge_max.unwrap())
                        {
                            for neighbor in graph.handle_edges_iter(handle, Direction::Right) {
                                let mut inc_kmer = kmer.clone();
                                inc_kmer.add_handle_to_complete(neighbor);

                                if next_count > 1 {
                                    inc_kmer.forks += 1;
                                }

                                incomplete_kmers.push(inc_kmer);
                            }
                        }
                    }
                }

                // Then complete all incomplete kmers
                while let Some(mut incomplete_kmer) = incomplete_kmers.pop() {
                    handle = incomplete_kmer.handles.pop().unwrap();
                    handle_seq = graph.sequence(handle).into_string_lossy();
                    handle_length = handle_seq.len() as u64;

                    let end = min(k - (incomplete_kmer.seq.len() as u64), handle_length);
                    let str_to_add = handle_seq.substring(0, end as usize).to_string();
                    incomplete_kmer.extend_kmer(str_to_add, handle);

                    if (incomplete_kmer.seq.len() as u64) == k {
                        if !complete_kmers.contains(&incomplete_kmer) {
                            complete_kmers.push(incomplete_kmer);
                        }
                    } else {
                        // NOTE: if there is no neighbor, the kmer does not get re-added
                        // to the incomplete ones, so that the external loop can end
                        for neighbor in graph.handle_edges_iter(handle, Direction::Right) {
                            let mut next_count: u64 = 0;
                            if edge_max.is_some() || degree_max.is_some() {
                                graph
                                    .handle_edges_iter(handle, Direction::Right)
                                    .for_each(|_| {
                                        next_count += 1;
                                    });
                            }

                            if (degree_max.is_none() && edge_max.is_none())
                                || (degree_max.is_some() && next_count < degree_max.unwrap())
                                || (edge_max.is_some() && incomplete_kmer.forks < edge_max.unwrap())
                            {
                                let mut inc_kmer = incomplete_kmer.clone();
                                inc_kmer.add_handle_to_complete(neighbor);

                                if next_count > 1 {
                                    inc_kmer.forks += 1;
                                }

                                incomplete_kmers.push(inc_kmer);
                            }
                        }
                    }
                }

                assert!(incomplete_kmers.is_empty());
            }
        }
    }

    complete_kmers
}

/// Struct that represents the kmer positions on the forward
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct KmerPos {
    //seq : String,
    /// The start position of the kmer
    pub(crate) start: u64,
    /// The end position of the kmer
    pub(crate) end: u64,
    /// The orientation of the kmer
    pub(crate) orient: bool,
}

fn get_seq_pos(
    handle: &Handle,
    node_ref: &Vec<NodeRef>,
    ref_length: &u64,
    handle_length: &u64,
) -> u64 {
    let pos: u64;

    let handle_rank: u64 = handle.unpack_number() - 1;
    let node_starting_pos: u64 = node_ref.get(handle_rank as usize).unwrap().seq_idx;

    if handle.is_reverse() {
        pos = ref_length - node_starting_pos - handle_length;
    } else {
        pos = node_starting_pos;
    }

    pos
}

// TODO: remove first parameter
pub fn generate_pos_on_ref(
    graph: &HashGraph,
    kmers_on_graph: &Vec<Kmer>,
    seq_length: &u64,
    node_ref: &Vec<NodeRef>,
) -> Vec<KmerPos> {
    let mut kmers_on_ref: Vec<KmerPos> = Vec::new();

    for kmer in kmers_on_graph {
        let first_handle_of_kmer = kmer.handles.first().unwrap();
        let first_handle_node = graph.get_node(&first_handle_of_kmer.id()).unwrap();
        let first_handle_length = first_handle_node.sequence.len();

        let last_handle_of_kmer = kmer.handles.last().unwrap();
        let last_handle_node = graph.get_node(&last_handle_of_kmer.id()).unwrap();
        let last_handle_length = last_handle_node.sequence.len();

        let start_ref = get_seq_pos(
            first_handle_of_kmer,
            node_ref,
            seq_length,
            &(first_handle_length as u64),
        ) + kmer.begin;
        let end_ref = get_seq_pos(
            last_handle_of_kmer,
            node_ref,
            seq_length,
            &(last_handle_length as u64),
        ) + kmer.end;

        let mut pos = KmerPos {
            start: start_ref,
            end: end_ref,
            orient: kmer.handle_orient,
        };
        kmers_on_ref.push(pos);
    }

    kmers_on_ref
}

/// This function converts kmer positions on the graph to kmer positions on the forward
/// (this is effectively a Kmer -> KmerPos conversion)
pub fn generate_pos_on_forward(
    kmers_on_graph: &Vec<Kmer>,
    forward: &String,
    seq_bw: &BitVec,
    node_ref: &Vec<NodeRef>,
) -> Vec<KmerPos> {
    let mut kmers_on_fwd: Vec<KmerPos> = Vec::new();

    for kmer in kmers_on_graph {
        // -1 required because i-eth node id is i-1-eth in node list

        let kmer_handle = kmer.handles.get(0).unwrap();
        let kmer_nodeId = kmer_handle.unpack_number() - 1;

        let noderef_kmer: &NodeRef = node_ref.get(kmer_nodeId as usize).unwrap();
        let start_pos_on_fwd = noderef_kmer.seq_idx + kmer.begin;
        let end_pos_on_fwd = noderef_kmer.seq_idx + kmer.end;

        let curr_kmer: KmerPos = KmerPos {
            //seq : kmer.seq.clone(),
            start: start_pos_on_fwd,
            end: end_pos_on_fwd,
            orient: true,
        };
        kmers_on_fwd.push(curr_kmer);
    }

    kmers_on_fwd
}

/// This function returns a HashMap representing the kmers and their positions,
/// having as key the kmer seq, and as values a KmerPos
pub fn generate_kmers_hash(
    kmers_on_graph: Vec<Kmer>,
    kmers_on_fwd: Vec<KmerPos>,
) -> AHashMap<String, KmerPos> {
    let mut kmers_hashed: AHashMap<String, KmerPos> = AHashMap::new();

    assert_eq!(kmers_on_graph.len(), kmers_on_fwd.len());

    for i in 0..kmers_on_graph.len() {
        let kmer_graph_ref = kmers_on_graph.get(i).unwrap();
        let seq = kmer_graph_ref.seq.clone();

        let pos_on_fwd = kmers_on_fwd.get(i).unwrap();

        //kmers_hashed.insert(seq, KmerPos { start: pos_on_fwd.start, end: pos_on_fwd.end});
    }

    kmers_hashed
}

pub fn generate_hash(kmers_on_graph: &Vec<Kmer>, hash_builder: &RandomState) -> Vec<u64> {
    let mut hashes: Vec<u64> = Vec::new();

    kmers_on_graph.iter().for_each(|kmer| {
        hashes.push(u64::get_hash(&kmer.seq, hash_builder));
    });

    hashes.sort();

    hashes
}
