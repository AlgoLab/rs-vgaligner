use handlegraph::hashgraph::HashGraph;
use handlegraph::handlegraph::HandleGraph;
use handlegraph::handle::{Handle, Direction, NodeId};
use std::cmp::min;
use boomphf::*;
use substring::Substring;
use ahash::AHashMap;
use crate::utils::NodeRef;
use bv::BitVec;
use std::collections::VecDeque;
use std::hash::BuildHasher;
use ahash::RandomState;
use ahash::CallHasher;

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
    pub(crate) handle: Handle
}

impl PartialEq for Kmer {
    fn eq(&self, other: &Self) -> bool {
        self.seq == other.seq
    }
}

impl Kmer {
    /// This function allows for a kmer to be extended, used when the
    /// kmer is on more than one handle
    fn extend_kmer(&mut self, new_seq : String) {
        self.seq.push_str(&new_seq);
        self.end += new_seq.len() as u64;
    }
}

/// This function computes a list of the unique kmers in a given HashGraph
pub fn generate_kmers(graph : &HashGraph, k : u64, edge_max : Option<u64>, degree_max : Option<u64>) -> Vec<Kmer> {

    let mut kmers : Vec<Kmer> = Vec::new();

    let mut prev_kmers_to_complete : VecDeque<Kmer> = VecDeque::new();

    // Sort both handles and edges since this is a PO
    let mut graph_handles : Vec<Handle> = graph.handles_iter().collect();
    graph_handles.sort();

    // For each handle
    graph_handles.iter().for_each(|h| {

        // Check both forward and reverse
        //for handle_is_rev in &[true, false] {

        //TODO: add reverse later
        for handle_is_rev in &[true] {
            let handle : Handle;

            let mut curr_kmers_to_complete : Vec<Kmer> = Vec::new();

            assert_eq!(curr_kmers_to_complete.len(), 0);

            match handle_is_rev {
                true => handle = *h,
                false => handle = h.flip()
            }

            if let Some(degree_max) = degree_max {
                let mut curr_count : u64 = 0;
                graph.handle_edges_iter(handle, Direction::Left).for_each(|_| curr_count += 1);
                if curr_count > degree_max {
                    continue;
                }
            }

            // Get current handle
            let node = graph.get_node(&handle.id()).unwrap();
            let handle_seq = node.sequence.to_string();
            let handle_length = handle_seq.len() as u64;

            // First try completing previously uncompleted kmers
            // EXAMPLE: suppose k=4
            // node i-1 = AC, node i = CGT
            // incomplete kmer is AC (len 2), I want ACCG (len 4 = k)
            // also C + CGT

            while !prev_kmers_to_complete.is_empty() {

                let mut incomplete_kmer = prev_kmers_to_complete.pop_front().unwrap();

                let end = min(k - (incomplete_kmer.seq.len() as u64), handle_length);
                let str_to_add = handle_seq.substring(0, end as usize).to_string();
                incomplete_kmer.extend_kmer(str_to_add);

                if (incomplete_kmer.seq.len() as u64) == k {

                    if !kmers.contains(&incomplete_kmer) {
                        kmers.push(incomplete_kmer);
                    }

                } else {
                    curr_kmers_to_complete.push(incomplete_kmer);
                }

            }

            // Then try to generate current node kmers

            if let Some(degree_max) = degree_max {
                let mut curr_count : u64 = 0;
                graph.handle_edges_iter(handle, Direction::Right).for_each(|neighbour|{
                    curr_count+=1;
                });
                if curr_count > degree_max {
                    continue;
                }
            }

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

            // Add kmers not completed in this iteration
            curr_kmers_to_complete.iter()
                .for_each(|inc_kmer| prev_kmers_to_complete.push_back(inc_kmer.clone()));
        }

    });

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

        let kmer_handle = kmer.handle;
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