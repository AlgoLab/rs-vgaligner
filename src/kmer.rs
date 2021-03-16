use handlegraph::hashgraph::HashGraph;
use handlegraph::handlegraph::HandleGraph;
use handlegraph::handle::{Handle, Direction};
use std::cmp::min;
use boomphf::*;
use substring::Substring;
use ahash::AHashMap;

#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub struct Kmer {
    seq : String,
    begin : u64,
    end : u64,
    handle: Handle
    //handle : Vec<Handle>,
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

    // For each handle
    graph.handles_iter().for_each(|h| {

        // Check both forward and reverse
        for handle_is_rev in &[true, false] {
            let handle : Handle;

            match handle_is_rev {
                true => handle = h,
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

            for i in 0..k {
                let begin = i;
                let end = min(i+k, handle_length);
                let mut kmer = Kmer {
                    seq : handle_seq.substring(begin as usize, end as usize).to_string(),
                    begin,
                    end,
                    handle
                    //handle : vec![handle]
                };

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

pub struct KmerPos {
    seq : String,
    start : u64,
    end : u64
}

pub fn generate_kmers_hash(kmers : &Vec<Kmer>) -> AHashMap<String, KmerPos> {
    let mut kmers_hashed: AHashMap<String, KmerPos> = AHashMap::new();

    for kmer in kmers {
        let pos = KmerPos {
            seq: kmer.seq.clone(),
            start: kmer.begin,
            end: kmer.end
        };
        kmers_hashed.insert(kmer.seq.clone(), pos);
    }

    

    kmers_hashed
}