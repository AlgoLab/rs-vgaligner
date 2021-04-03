use handlegraph::hashgraph::HashGraph;
use handlegraph::handlegraph::HandleGraph;
use handlegraph::handle::{Handle, Direction, NodeId};
use std::cmp::min;
use boomphf::*;
use substring::Substring;
use ahash::AHashMap;
use crate::utils::NodeRef;
use bv::BitVec;
use std::collections::{VecDeque, HashMap};
use std::hash::BuildHasher;
use ahash::RandomState;
use ahash::CallHasher;
use bstr::ByteVec;

/// Struct that represents a kmer in the graph
#[derive(Debug, Clone)]
pub struct Kmer {
    /// The sequence of the kmer
    pub(crate) seq : String,
    /// The start position relative to the handle
    begin : u64,
    /// The end position relative to the handle
    end : u64,
    /// The handle where the kmer starts (NOTE: this is only the start
    /// the actual kmer may be on more than on handle)
    pub(crate) handles : Vec<Handle>
}

impl PartialEq for Kmer {
    fn eq(&self, other: &Self) -> bool {
        self.seq == other.seq
    }
}

impl Kmer {
    /// This function allows for a kmer to be extended, used when the
    /// kmer is on more than one handle
    fn extend_kmer(&mut self, new_seq : String, new_handle : Handle) {
        self.seq.push_str(&new_seq);
        self.end += new_seq.len() as u64;
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

pub fn generate_kmers(graph : &HashGraph, k : u64, edge_max : Option<u64>, degree_max : Option<u64>) -> Vec<Kmer> {
    let mut complete_kmers : Vec<Kmer> = Vec::new();

    // Sort both handles
    let mut sorted_graph_handles : Vec<Handle> = graph.handles_iter().collect();
    sorted_graph_handles.sort();

    let mut incomplete_kmers : HashMap<Handle, VecDeque<Kmer>> = HashMap::with_capacity(sorted_graph_handles.len());
    sorted_graph_handles.iter().for_each(|h| { incomplete_kmers.insert(h.clone(), VecDeque::new()); });

    for handle in sorted_graph_handles {

        // Get current node/handle
        let node = graph.get_node(&handle.id()).unwrap();
        let handle_seq = node.sequence.to_string();
        let handle_length = handle_seq.len() as u64;

        // First try completing incomplete kmers
        let mut curr_incomplete_kmers = incomplete_kmers.get_mut(&handle).unwrap().clone();
        while let Some(mut incomplete_kmer) = curr_incomplete_kmers.pop_front() {

            let end = min(k - (incomplete_kmer.seq.len() as u64), handle_length);
            let str_to_add = handle_seq.substring(0, end as usize).to_string();
            incomplete_kmer.extend_kmer(str_to_add, handle);

            if (incomplete_kmer.seq.len() as u64) == k {
                if !complete_kmers.contains(&incomplete_kmer) {
                    complete_kmers.push(incomplete_kmer);
                }

            } else {
                for neighbor in graph.handle_edges_iter(handle, Direction::Right) {
                    incomplete_kmers.get_mut(&neighbor).unwrap().push_back(incomplete_kmer.clone());
                }

            }
        }

        assert!(curr_incomplete_kmers.is_empty());
        incomplete_kmers.get_mut(&handle).unwrap().clear();

        // Then try generating the kmers from the given node
        for i in 0..handle_length {
            let begin = i;
            let end = min(i+k, handle_length);
            let mut kmer = Kmer {
                seq : handle_seq.substring(begin as usize, end as usize).to_string(),
                begin,
                end,
                handles : vec![handle]
            };

            if (kmer.seq.len() as u64) == k {
                if !complete_kmers.contains(&kmer) {
                    complete_kmers.push(kmer);
                }
            } else {
                for neighbor in graph.handle_edges_iter(handle, Direction::Right) {
                    incomplete_kmers.get_mut(&neighbor).unwrap().push_back(kmer.clone());
                }
            }
        }

    }

    complete_kmers
}

pub fn generate_kmers_rev(graph : &HashGraph, k : u64, edge_max : Option<u64>, degree_max : Option<u64>) -> Vec<Kmer> {
    let mut complete_kmers : Vec<Kmer> = Vec::new();

    // Sort both handles
    let mut sorted_graph_handles : Vec<Handle> = graph.handles_iter().collect();
    sorted_graph_handles.sort();
    sorted_graph_handles.reverse();

    let mut incomplete_kmers : HashMap<Handle, VecDeque<Kmer>> = HashMap::with_capacity(sorted_graph_handles.len());
    sorted_graph_handles.iter().for_each(|h| { incomplete_kmers.insert(h.flip().clone(), VecDeque::new()); });

    for handle_forward in sorted_graph_handles {

        let handle = handle_forward.flip();

        // Get current node/handle
        let handle_seq = graph.sequence(handle).into_string_lossy();
        let handle_length = handle_seq.len() as u64;

        // First try completing incomplete kmers
        let mut curr_incomplete_kmers = incomplete_kmers.get_mut(&handle).unwrap().clone();
        while let Some(mut incomplete_kmer) = curr_incomplete_kmers.pop_front() {

            let end = min(k - (incomplete_kmer.seq.len() as u64), handle_length);
            let str_to_add = handle_seq.substring(0, end as usize).to_string();
            incomplete_kmer.extend_kmer(str_to_add, handle);

            if (incomplete_kmer.seq.len() as u64) == k {
                if !complete_kmers.contains(&incomplete_kmer) {
                    complete_kmers.push(incomplete_kmer);
                }

            } else {
                // TODO: kmer should contain multiple handles
                for neighbor in graph.handle_edges_iter(handle, Direction::Right) {
                    incomplete_kmers.get_mut(&neighbor).unwrap().push_back(incomplete_kmer.clone());
                }

            }
        }

        assert!(curr_incomplete_kmers.is_empty());
        incomplete_kmers.get_mut(&handle).unwrap().clear();

        // Then try generating the kmers from the given node
        for i in 0..handle_length {
            let begin = i;
            let end = min(i+k, handle_length);
            let mut kmer = Kmer {
                seq : handle_seq.substring(begin as usize, end as usize).to_string(),
                begin,
                end,
                handles : vec![handle]
            };

            if (kmer.seq.len() as u64) == k {
                if !complete_kmers.contains(&kmer) {
                    complete_kmers.push(kmer);
                }
            } else {
                // TODO: kmer should contain multiple handles
                for neighbor in graph.handle_edges_iter(handle, Direction::Right) {
                    incomplete_kmers.get_mut(&neighbor).unwrap().push_back(kmer.clone());
                }
            }
        }

    }

    complete_kmers
}

pub fn merge_kmers(kmers_fwd : Vec<Kmer>, kmers_rev : Vec<Kmer>) -> Vec<Kmer> {
    let mut kmers : Vec<Kmer> = Vec::new();

    kmers = kmers_fwd.clone();

    for kmer in kmers_rev {
        if !kmers.contains(&kmer) {
            kmers.push(kmer);
        }
    }

    kmers
}

/// Struct that represents the kmer positions on the forward
#[derive(Debug)]
pub struct KmerPos {
    //seq : String,
    /// The start position of the kmer on the forward
    pub(crate) start : u64,
    /// The end position of the kmer on the forward
    pub(crate) end : u64
}

/// This function converts kmer positions on the graph to kmer positions on the forward
/// (this is effectively a Kmer -> KmerPos conversion)
pub fn generate_pos_on_forward(kmers_on_graph : &Vec<Kmer>, forward : &String,
                               seq_bw : &BitVec, node_ref : &Vec<NodeRef>) -> Vec<KmerPos> {
    let mut kmers_on_fwd : Vec<KmerPos> = Vec::new();

    for kmer in kmers_on_graph {

        // -1 required because i-eth node id is i-1-eth in node list
        // see TODO in kmer

        let kmer_handle = kmer.handles.get(0).unwrap();
        let kmer_nodeId = kmer_handle.unpack_number()-1;

        let noderef_kmer : &NodeRef = node_ref.get(kmer_nodeId as usize).unwrap();
        let start_pos_on_fwd = noderef_kmer.seq_idx + kmer.begin;
        let end_pos_on_fwd = noderef_kmer.seq_idx + kmer.end;

        let curr_kmer : KmerPos = KmerPos {
            //seq : kmer.seq.clone(),
            start : start_pos_on_fwd,
            end : end_pos_on_fwd
        };
        kmers_on_fwd.push(curr_kmer);
    }

    kmers_on_fwd
}

/// This function returns a HashMap representing the kmers and their positions,
/// having as key the kmer seq, and as values a KmerPos
pub fn generate_kmers_hash(kmers_on_graph : Vec<Kmer>, kmers_on_fwd : Vec<KmerPos>) -> AHashMap<String, KmerPos> {
    let mut kmers_hashed: AHashMap<String, KmerPos> = AHashMap::new();

    assert_eq!(kmers_on_graph.len(), kmers_on_fwd.len());

    for i in 0..kmers_on_graph.len() {

        let kmer_graph_ref = kmers_on_graph.get(i).unwrap();
        let seq = kmer_graph_ref.seq.clone();

        let pos_on_fwd  = kmers_on_fwd.get(i).unwrap();

        //kmers_hashed.insert(seq, KmerPos { start: pos_on_fwd.start, end: pos_on_fwd.end});
    }
    

    kmers_hashed
}

pub fn generate_hash(kmers_on_graph : &Vec<Kmer>, hash_builder : &RandomState) -> Vec<u64> {
    let mut hashes : Vec<u64> = Vec::new();

    kmers_on_graph.iter().for_each(|kmer|{
        hashes.push(u32::get_hash(&kmer.seq, hash_builder));
    });

    hashes
}

pub fn generate_mphf(kmer_hashes : &Vec<u64>) -> Mphf<u64> {
    let phf = Mphf::new(1.7, &kmer_hashes.clone());
    phf
}