use handlegraph::hashgraph::HashGraph;
use handlegraph::handlegraph::HandleGraph;
use handlegraph::handle::{Handle, Direction};
use std::cmp::min;
use boomphf::*;
use substring::Substring;

#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub struct kmer {
    seq : String,
    begin : u64,
    end : u64,
    handle : Handle,
}

impl kmer {
    pub fn extend_kmer(&mut self, new_seq : String) {
        self.seq.push_str(&new_seq);
        self.end += new_seq.len() as u64;
        // This makes the end go outside the handle...
        // TODO: rethink how kmers are stored (maybe Vec<Handle> and multiple starts/ends?)
    }
}

// TODO: currently only returns kmers in the same node
pub fn generate_kmers(graph : &HashGraph, k : u64, degree_max : Option<u64>) -> Vec<kmer> {

    let mut kmers : Vec<kmer> = Vec::new();

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
                let kmer = kmer {
                    // This cannot be right...
                    //seq : (&handle_seq.to_owned()[begin as usize..end as usize]).to_string(),
                    // This is better, but requires an external crate
                    seq : handle_seq.substring(begin as usize, end as usize).to_string(),
                    begin,
                    end,
                    handle
                };

                if kmer.seq.len() < k as usize {
                    // TODO: check neighbors
                    // Will do it later, once I understand how kmers are used
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

pub fn generate_kmers_hash(kmers : &Vec<kmer>) -> Vec<u64> {
    let phf = Mphf::new(1.7, kmers.clone().as_slice());

    // Get hash value of all objects
    let mut hashes = Vec::new();
    for v in kmers {
        hashes.push(phf.hash(&v));
    }
    hashes.sort();

    hashes
}