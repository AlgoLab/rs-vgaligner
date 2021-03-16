use handlegraph::hashgraph::HashGraph;
use handlegraph::handlegraph::HandleGraph;
use handlegraph::handle::{Handle, Direction, NodeId};
use std::cmp::min;
use boomphf::*;
use substring::Substring;
use ahash::AHashMap;
use crate::utils::NodeRef;
use bv::BitVec;

#[derive(Debug)]
pub struct Kmer {
    seq : String,
    begin : u64,
    end : u64,
    handle: Handle
    //handle : Vec<Handle>,
}

impl PartialEq for Kmer {
    fn eq(&self, other: &Self) -> bool {
        self.seq == other.seq
    }
}

/*
impl Kmer {
    pub fn extend_kmer(&mut self, new_seq : String, new_handle : Handle) {
        self.seq.push_str(&new_seq);
        self.end += new_seq.len() as u64;
        self.handle.push(new_handle);
    }
}
*/

pub fn generate_kmers(graph : &HashGraph, k : u64, degree_max : Option<u64>) -> Vec<Kmer> {

    let mut kmers : Vec<Kmer> = Vec::new();

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

            match handle_is_rev {
                true => handle = *h,
                false => handle = h.flip()
            }

            if let Some(degree_max) = degree_max {
                let mut curr_count : u64 = 0;
                graph.handle_edges_iter(handle, Direction::Right).for_each(|neighbour|{
                    curr_count+=1;
                });
                if curr_count > degree_max {
                    continue;
                }
            }

            let node = graph.get_node(&handle.id()).unwrap();
            let handle_seq = node.sequence.to_string();
            let handle_length = handle_seq.len() as u64;

            // Can't generate single-node kmer if the seq is shorter than k
            if handle_length < k {
                continue
            }

            // NOTE: end is non-inclusive
            for i in 0..(handle_length-k+1) {
                let begin = i;
                let end = min(i+k, handle_length);
                let mut kmer = Kmer {
                    seq : handle_seq.substring(begin as usize, end as usize).to_string(),
                    begin,
                    end,
                    handle
                    //handle : vec![handle]
                };
                //println!("Seq: {} Begin: {}, End: {}, Sub: {}", handle_seq,begin,end, kmer.seq);

                if kmer.seq.len() < k as usize {
                    /*
                    // maybe too taxing...
                    // also requires some kind of graph visit to be 100% correct
                    for neighbour in graph.handle_edges_iter(handle, Direction::Right) {

                        let neighbor_node = graph.get_node(&neighbour.id()).unwrap();
                        let neighbor_seq = neighbor_node.sequence.to_string();
                        let neighbor_length = neighbor_seq.len() as u64;

                        let remaining_len = k - kmer.seq.len() as u64;
                        let extension_len = min(remaining_len, neighbor_length);

                        let mut ext_kmer = kmer.clone();

                        let seq_to_add = neighbor_seq.substring(0 as usize, extension_len as usize).to_string();

                        ext_kmer.extend_kmer(seq_to_add, neighbour);

                        if ext_kmer.seq.len() == k as usize && !kmers.contains(&ext_kmer) {
                            kmers.push(ext_kmer);
                        }

                    }
                    */
                } else {
                    // Kmers must be unique for hashing
                    if !kmers.contains(&kmer) {
                        kmers.push(kmer);
                    }
                }
            }

        }

    });

    kmers
}

/* TODO
pub fn extend_kmers_with_BFS(graph : &HashGraph, k : u64, startHandle : &Handle, curr_kmer : &mut Kmer) ->  {
    let mut extended_kmers : Vec<Kmer> = Vec::new();
    extended_kmers
}
 */

#[derive(Debug)]
pub struct KmerPos {
    //seq : String,
    start : u64,
    end : u64
}

pub fn generate_pos_on_forward(kmers_on_graph : &Vec<Kmer>, forward : &String,
                               seq_bw : &BitVec, node_ref : &Vec<NodeRef>  ) -> Vec<KmerPos> {
    let mut kmers_on_fwd : Vec<KmerPos> = Vec::new();

    //println!("Node_ref: {:#?}", node_ref);
    //println!("Len: {}", node_ref.len());
    for kmer in kmers_on_graph {

        // -1 required because i-eth node id is i-1-eth in node list
        let kmer_nodeId = kmer.handle.unpack_number()-1;
        //println!("id: {} value: {}",kmer_nodeId, kmer.seq);
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

    //kmers_on_fwd.sort_by(|a, b| a.start.cmp(&b.start));

    kmers_on_fwd
}

pub fn generate_kmers_hash(kmers_on_graph : Vec<Kmer>, kmers_on_fwd : Vec<KmerPos>) -> AHashMap<String, KmerPos> {
    let mut kmers_hashed: AHashMap<String, KmerPos> = AHashMap::new();

    for i in 0..kmers_on_graph.len() {

        let kmer_graph_ref = kmers_on_graph.get(i).unwrap();
        let seq = kmer_graph_ref.seq.clone();

        let pos_on_fwd  = kmers_on_fwd.get(i).unwrap();

        kmers_hashed.insert(seq, KmerPos { start: pos_on_fwd.start, end: pos_on_fwd.end});
    }
    

    kmers_hashed
}