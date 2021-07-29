use std::cmp::{min, Ordering};
use std::collections::VecDeque;

use ahash::CallHasher;
use ahash::RandomState;

use bstr::ByteVec;
use handlegraph::handle::{Direction, Handle};
use handlegraph::handlegraph::HandleGraph;
use handlegraph::hashgraph::HashGraph;
use handlegraph::pathgraph::PathHandleGraph;
use serde::{Deserialize, Serialize};
use substring::Substring;

use crate::serialization::SerializableHandle;
use crate::utils::NodeRef;
use std::ops::{Index, Deref};
use rayon::prelude::ParallelSliceMut;

/// Represents a kmer in the graph
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Kmer {
    /// The sequence of the kmer
    pub(crate) seq: String,
    /// The start position relative to the handle
    pub(crate) begin: u64,
    /// The end position relative to the handle
    pub(crate) end: u64,
    /// The first handle of the kmer
    #[serde(with = "SerializableHandle")]
    pub(crate) first: Handle,
    /// The last handle of the kmer
    #[serde(with = "SerializableHandle")]
    pub(crate) last: Handle,
    /// The orientation of the handles
    pub(crate) handle_orient: bool,
    /// The number of forks the kmer has been involved in
    pub(crate) forks: u64,
}

/*
impl PartialEq for Kmer {
    fn eq(&self, other: &Self) -> bool {
        self.seq == other.seq && self.begin == other.begin && self.end == other.end
        && self.first == other.first && self.last == other.last &&
        self.handle_orient == other.handle_orient
    }
}
 */

impl Kmer {
    /// This function allows for a kmer to be extended, used when the
    /// kmer is on more than one handle
    fn extend_kmer(&mut self, new_seq: String, new_handle: Handle) {
        self.seq.push_str(&new_seq);
        //self.end += new_seq.len() as u64;
        self.end = new_seq.len() as u64;
        self.last = new_handle;
    }

    fn add_handle_to_complete(&mut self, new_handle: Handle) {
        self.last = new_handle;
    }
}

/// Generate kmers having size k from a given HashGraph. Both the number of overall visited edges
/// and maximum (outgoing) degree for a node can be limited, which should help with bigger graphs.
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

        // Kmers will be generated from both forward and reverse, as handlegraph is capable
        // of storing both
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

            // Check if the handle has more outgoing edges than the maximum allowed,
            // and if that's the case, skip the current handle
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

            // This will store kmers that have not yet reached size k, which
            // will have to be completed from the neighbours of the handle
            let mut incomplete_kmers: Vec<Kmer> = Vec::new();

            // Try generating the "internal" kmers from the given handle
            for i in 0..handle_length {
                let begin = i;
                let end = min(i + k, handle_length);
                let kmer = Kmer {
                    seq: handle_seq
                        .substring(begin as usize, end as usize)
                        .to_string(),
                    begin,
                    end,
                    first: handle,
                    last: handle,
                    handle_orient: orient,
                    forks: 0,
                };

                // Ignore Ns in kmer generation
                if kmer.seq.contains('N') {
                    continue;
                }

                // If the kmer has already reached size k...
                // NOTE: this implies that the sequence encoded by the current handle
                // has size >= k
                if (kmer.seq.len() as u64) == k {
                    //if !complete_kmers.contains(&kmer) {
                        complete_kmers.push(kmer);
                    //}
                } else {
                    // The kmer is incomplete, thus will have to be completed to reach size k

                    // Check that eventual limits have not been reach yet, otherwise ignore
                    // this kmer

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

                        // Create a copy of the incomplete kmer for each neighbour handle,
                        // so that they can be completed

                        for neighbor in graph.handle_edges_iter(handle, Direction::Right) {
                            let mut inc_kmer = kmer.clone();
                            inc_kmer.last = neighbor;

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
                handle = incomplete_kmer.last;
                handle_seq = graph.sequence(handle).into_string_lossy();
                handle_length = handle_seq.len() as u64;

                let end = min(k - (incomplete_kmer.seq.len() as u64), handle_length);
                let str_to_add = handle_seq.substring(0, end as usize).to_string();
                incomplete_kmer.extend_kmer(str_to_add, handle);

                // Ignore Ns during kmer generation
                if incomplete_kmer.seq.contains('N') {
                    continue;
                }

                if (incomplete_kmer.seq.len() as u64) == k {
                    //if !complete_kmers.contains(&incomplete_kmer) {
                        complete_kmers.push(incomplete_kmer);
                    //}
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

            // IMPORTANT NOTE: a single iteration (of the most external for loop) completes
            // ALL possible kmers that start in the given handle. If the graph "ends" but the
            // kmer cannot reach size k (i.e. it is incomplete) it gets discarded.
            assert!(incomplete_kmers.is_empty());
        }
    }

    // Sort the kmers so that equal kmers (= having the same sequence) are close to each other
    // Note that the same kmer can appear in different places
    // (e.g. CACTTCAC -> CAC and CAC must be consecutive in the ordering)
    complete_kmers.par_sort_by(|x,y| x.seq.cmp(&y.seq));
    // Also dedup the vec as exact duplicates only waste space. Also note that dedup only works
    // on consecutive duplicates, so only by sorting beforehand it works correctly.
    complete_kmers.dedup();

    complete_kmers
}

/// Generate kmers having size k from a given HashGraph, using paths to guide the generation. This
/// should require less memory than generate_kmers, but may not return the same kmers (it depends
/// on how paths are defined in the GFA).
pub fn generate_kmers_linearly(
    graph: &HashGraph,
    k: u64,
    edge_max: Option<u64>,
    degree_max: Option<u64>,
) -> Vec<Kmer> {
    assert!(!graph.paths.is_empty());

    // This requires two iterations, one on the forward (handles) and one on the reverse (handles).
    let mut forward_kmers = generate_kmers_linearly_forward(graph, k, edge_max, degree_max);
    let mut reverse_kmers = generate_kmers_linearly_reverse(graph, k, edge_max, degree_max);

    // Merge the kmers obtained previously
    let mut kmers = merge_kmers(forward_kmers, reverse_kmers);

    // Sort the kmers so that equal kmers (= having the same sequence) are close to each other
    // Note that the same kmer can appear in different places
    // (e.g. CACTTCAC -> CAC and CAC must be consecutive in the ordering)
    kmers.par_sort_by(|x,y| x.seq.cmp(&y.seq));
    // Also dedup the vec as exact duplicates only waste space. Also note that dedup only works
    // on consecutive duplicates, so only by sorting beforehand it works correctly.
    kmers.dedup();

    kmers
}

/// Generate kmers having size k from the forward strand encoded in a HashGraph
fn generate_kmers_linearly_forward(
    graph: &HashGraph,
    k: u64,
    _edge_max: Option<u64>,
    _degree_max: Option<u64>,
) -> Vec<Kmer> {
    let mut kmers: Vec<Kmer> = Vec::new();
    let mut prev_kmers_to_complete: VecDeque<Kmer> = VecDeque::new();

    for path_id in graph.paths_iter() {
        let path = graph.get_path(path_id).unwrap();

        //println!("Path: {:#?}", path);
        prev_kmers_to_complete.clear();

        // Navigate each path linearly, i.e. in the order of the nodes that make up the path
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

                // Ignore Ns during kmer generation
                if incomplete_kmer.seq.contains('N') {
                    continue;
                }

                if (incomplete_kmer.seq.len() as u64) == k {
                    //if !kmers.contains(&incomplete_kmer) {
                        kmers.push(incomplete_kmer);
                    //}
                } else {
                    curr_kmers_to_complete.push(incomplete_kmer);
                }
            }

            // Then try to generate current node kmers
            for i in 0..handle_length {
                let begin = i;
                let end = min(i + k, handle_length);
                let kmer = Kmer {
                    seq: handle_seq
                        .substring(begin as usize, end as usize)
                        .to_string(),
                    begin,
                    end,
                    first: *handle,
                    last: *handle,
                    handle_orient: true,
                    forks: 0,
                };

                // Ignore Ns during kmer generation
                if kmer.seq.contains('N') {
                    continue;
                }

                if (kmer.seq.len() as u64) == k {
                    //if !kmer.seq.contains('N') && !kmers.contains(&kmer) {
                        kmers.push(kmer);
                    //}
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

/// Generate kmers having size k from the reverse strand encoded in a HashGraph
pub fn generate_kmers_linearly_reverse(
    graph: &HashGraph,
    k: u64,
    _edge_max: Option<u64>,
    _degree_max: Option<u64>,
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

                // Ignore Ns during kmer generation
                if incomplete_kmer.seq.contains('N') {
                    continue;
                }

                if (incomplete_kmer.seq.len() as u64) == k {
                    //if !kmers.contains(&incomplete_kmer) {
                        kmers.push(incomplete_kmer);
                    //}
                } else {
                    curr_kmers_to_complete.push(incomplete_kmer);
                }
            }

            // Then try to generate current node kmers
            for i in 0..handle_length {
                let begin = i;
                let end = min(i + k, handle_length);
                let kmer = Kmer {
                    seq: handle_seq
                        .substring(begin as usize, end as usize)
                        .to_string(),
                    begin,
                    end,
                    first: handle,
                    last: handle,
                    handle_orient: false,
                    forks: 0,
                };

                // Ignore Ns during kmer generation
                if kmer.seq.contains('N') {
                    continue;
                }

                if (kmer.seq.len() as u64) == k {
                    //if !kmers.contains(&kmer) {
                        kmers.push(kmer);
                    //}
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

/// Merge the kmers obtained by generate_kmers_linearly_forward and generate_kmers_linearly_reverse
fn merge_kmers(kmers_fwd: Vec<Kmer>, kmers_rev: Vec<Kmer>) -> Vec<Kmer> {
    let mut kmers: Vec<Kmer> = kmers_fwd.clone();

    for kmer in kmers_rev {
        //if !kmers.contains(&kmer) {
            kmers.push(kmer);
        //}
    }

    kmers
}

/// Represent kmer positions on the forward/reverse linearization
#[derive(Debug, Clone, Eq, PartialEq, Serialize, Deserialize)]
pub struct KmerPos {
    //seq : String,
    /// The start position of the kmer
    pub(crate) start: u64,
    /// The end position of the kmer
    pub(crate) end: u64,
    /// The orientation of the kmer
    pub(crate) orient: bool,
}

/// Obtain the position of a given handle in the serialized forward/reverse
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

/// This function converts kmer positions on the graph to kmer positions on the forward/reverse
/// linearization (this is effectively a Kmer -> KmerPos conversion)
pub fn generate_pos_on_ref(
    graph: &HashGraph,
    kmers_on_graph: &Vec<Kmer>,
    seq_length: &u64,
    node_ref: &Vec<NodeRef>,
) -> Vec<KmerPos> {
    let mut kmers_on_ref: Vec<KmerPos> = Vec::new();

    for kmer in kmers_on_graph {
        let first_handle_of_kmer = kmer.first;
        let first_handle_node = graph.get_node(&first_handle_of_kmer.id()).unwrap();
        let first_handle_length = first_handle_node.sequence.len();

        let last_handle_of_kmer = kmer.last;
        let last_handle_node = graph.get_node(&last_handle_of_kmer.id()).unwrap();
        let last_handle_length = last_handle_node.sequence.len();

        let start_ref = get_seq_pos(
            &first_handle_of_kmer,
            node_ref,
            seq_length,
            &(first_handle_length as u64),
        ) + kmer.begin;
        let end_ref = get_seq_pos(
            &last_handle_of_kmer,
            node_ref,
            seq_length,
            &(last_handle_length as u64),
        ) + kmer.end;

        let pos = KmerPos {
            start: start_ref,
            end: end_ref,
            orient: kmer.handle_orient,
        };
        kmers_on_ref.push(pos);
    }

    kmers_on_ref
}

/// This function converts kmer positions on the graph to kmer positions on the forward/reverse
/// linearization (this is effectively a Kmer -> KmerPos conversion)
pub fn generate_pos_on_ref_2(
    graph: &HashGraph,
    kmers_on_graph: &Vec<Kmer>,
    seq_length: &u64,
    node_ref: &Vec<NodeRef>,
    hashes: &mut Vec<u64>,
    kmers_start_offsets: &mut Vec<u64>,
) -> Vec<KmerPos> {
    let mut kmers_on_ref: Vec<Vec<KmerPos>> = Vec::new();

    // Hash builder to generate hashes
    //let hash_builder = RandomState::with_seeds(0, 0, 0, 0);
    let mut last_kmer: Option<String> = None;
    let mut curr_kmer_positions : Vec<KmerPos> = Vec::new();

    for kmer in kmers_on_graph {

        let first_handle_of_kmer = kmer.first;
        let first_handle_node = graph.get_node(&first_handle_of_kmer.id()).unwrap();
        let first_handle_length = first_handle_node.sequence.len();

        let last_handle_of_kmer = kmer.last;
        let last_handle_node = graph.get_node(&last_handle_of_kmer.id()).unwrap();
        let last_handle_length = last_handle_node.sequence.len();

        let start_ref = get_seq_pos(
            &first_handle_of_kmer,
            node_ref,
            seq_length,
            &(first_handle_length as u64),
        ) + kmer.begin;
        let end_ref = get_seq_pos(
            &last_handle_of_kmer,
            node_ref,
            seq_length,
            &(last_handle_length as u64),
        ) + kmer.end;

        let pos = KmerPos {
            start: start_ref,
            end: end_ref,
            orient: kmer.handle_orient,
        };

        last_kmer = match last_kmer {
            // This represents the first iteration
            None => {
                curr_kmer_positions.push(pos);
                Some(kmer.clone().seq)
            },
            Some(last_kmer) => {
                // If the kmer has changed, push the hash of the previous one and the
                // list of positions
                if last_kmer != kmer.seq {
                    let kmer_hash  = generate_hash(&last_kmer);
                    hashes.push(kmer_hash);

                    // Push the positions of the current kmer
                    kmers_on_ref.push(curr_kmer_positions.clone());

                    // Start with new value
                    curr_kmer_positions = Vec::new();
                    curr_kmer_positions.push(pos);
                    Some(kmer.clone().seq)
                } else {
                    // if the kmer has not changed, just add it to the list of positions
                    curr_kmer_positions.push(pos);
                    Some(last_kmer)
                }
            }
        };

        if kmer == kmers_on_graph.last().unwrap() {
            let kmer_hash  = generate_hash(&kmer.seq);
            hashes.push(kmer_hash);
            kmers_on_ref.push(curr_kmer_positions.clone());
        }

    }

    // sort each kmer's positions
    // (dedup not necessary since already done when generating graph kmers)
    for positions_list in &mut kmers_on_ref {
        //println!("Unsorted positions: {:#?}", positions_list);
        //positions_list.dedup();
        positions_list.par_sort_by(|x,y| match x.orient.cmp(&y.orient) {
            Ordering::Equal => x.start.cmp(&y.start),
            other => other,
        });
        //println!("Sorted positions: {:#?}", positions_list);
    }
    //println!("Kmers on ref non flattened: {:#?}", kmers_on_ref);

    // flatten kmer_pos_on_ref and store offsets.
    // Instead of having separate vecs for each kmers
    // (e.g. [kmer1 : [pos11, pos12, ... pos1n], kmer2 : [pos21, pos22, ... pos2m] ... ])
    // obtain a single list, but also store where each kmer starts
    // (e.g. [pos11, pos12, ... pos1n, pos21, pos22, ... pos2m ...] <- kmers positions
    // and [0, n ...] <- offset where each kmer start)
    // this is coherent with the original C++ implementation

    let mut kmers_on_ref_flattened : Vec<KmerPos> = Vec::new();
    let mut offset : u64 = 0;

    for positions_list in kmers_on_ref {
        kmers_start_offsets.push(offset);
        for position in positions_list {
            kmers_on_ref_flattened.push(position);
            offset += 1;
        }

        // Add end delimiter, this is needed because
        // we only store the starting positions as offsets
        let delimiter = KmerPos {
            start: u64::max_value(),
            end: u64::max_value(),
            orient: false  //Doesn't really matter, will only check max_value
        };
        kmers_on_ref_flattened.push(delimiter);
        offset += 1;
    }
    //println!("Kmers on ref flattened: {:#?}", kmers_on_ref_flattened);
    //println!("Kmers start offset: {:#?}", kmers_start_offsets);

    kmers_on_ref_flattened
}

/*
// Not used anymore, may still be useful
fn generate_pos_on_forward(
    kmers_on_graph: &Vec<Kmer>,
    _forward: &String,
    _seq_bw: &BitVec,
    node_ref: &Vec<NodeRef>,
) -> Vec<KmerPos> {
    let mut kmers_on_fwd: Vec<KmerPos> = Vec::new();

    for kmer in kmers_on_graph {
        // -1 required because i-eth node id is i-1-eth in node list
        let kmer_handle = kmer.first;
        let kmer_node_id = kmer_handle.unpack_number() - 1;

        let noderef_kmer: &NodeRef = node_ref.get(kmer_node_id as usize).unwrap();
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
 */


/// Generate the hashes of the sequence encoded in each kmer
pub fn generate_hash(seq: &String) -> u64 {
    let hash_builder = RandomState::with_seeds(0, 0, 0, 0);
    u64::get_hash(&seq, &hash_builder)
}
