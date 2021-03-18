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

// TODO: I only need the first/starting handle, so a single handle
// would probably be better..
#[derive(Debug, Clone)]
pub struct Kmer {
    pub(crate) seq : String,
    begin : u64,
    end : u64,
    //handle: Handle
    handle : Vec<Handle>,
}

impl PartialEq for Kmer {
    fn eq(&self, other: &Self) -> bool {
        self.seq == other.seq
    }
}

impl Kmer {
    pub fn extend_kmer(&mut self, new_seq : String, new_handle : Handle) {
        self.seq.push_str(&new_seq);
        self.end += new_seq.len() as u64;
        self.handle.push(new_handle);
    }
}

pub fn generate_kmers(graph : &HashGraph, k : u64, degree_max : Option<u64>) -> Vec<Kmer> {

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
                incomplete_kmer.extend_kmer(str_to_add, handle);

                if (incomplete_kmer.seq.len() as u64) == k {

                    //if !kmers.contains(&incomplete_kmer) {
                        kmers.push(incomplete_kmer);
                    //}

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
                    //handle
                    handle : vec![handle]
                };

                if (kmer.seq.len() as u64) == k {
                    //if !kmers.contains(&kmer) {
                        kmers.push(kmer);
                    //}
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

#[derive(Debug)]
pub struct KmerPos {
    //seq : String,
    start : u64,
    end : u64
}

pub fn generate_pos_on_forward(kmers_on_graph : &Vec<Kmer>, forward : &String,
                               seq_bw : &BitVec, node_ref : &Vec<NodeRef>  ) -> Vec<KmerPos> {
    let mut kmers_on_fwd : Vec<KmerPos> = Vec::new();

    for kmer in kmers_on_graph {

        // -1 required because i-eth node id is i-1-eth in node list
        // see TODO in kmer

        let kmer_handle = kmer.handle.get(0).unwrap();
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