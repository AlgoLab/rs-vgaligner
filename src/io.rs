use crate::kmer::{Kmer, KmerPos};
use bv::BitVec;

/// This function prints kmers data
pub fn print_kmers(kmers_on_graph : &Vec<Kmer>, kmers_on_seq_fwd : &Vec<KmerPos>) {
    for i in 0..kmers_on_graph.len() {
        let graph_kmer = kmers_on_graph.get(i).unwrap();
        let fwd_kmer = kmers_on_seq_fwd.get(i).unwrap();

        println!("Seq: {},\nStart_fwd: {},\nEnd_fwd: {},\nHandles: {:#?}\n",
                 graph_kmer.seq, fwd_kmer.start, fwd_kmer.end, graph_kmer.handle);
    }
}

/// This function prints the bitvector
pub fn print_bitvec(bv : &BitVec) {
    for i in 0..bv.len() {
        let bit = bv.get(i);
        match bit {
            true => print!("{}",1),
            false => print!("{}",0),
        }
    }
}