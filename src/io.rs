use crate::kmer::{Kmer, KmerPos};
use bv::BitVec;
use std::path::PathBuf;
use std::fs::File;
use std::io::Write;
use std::ptr::null;

/// This function prints kmers data
pub fn print_kmers(kmers_on_graph : &Vec<Kmer>, kmers_on_seq_fwd : &Vec<KmerPos>) {
    for i in 0..kmers_on_graph.len() {
        let graph_kmer = kmers_on_graph.get(i).unwrap();
        let fwd_kmer = kmers_on_seq_fwd.get(i).unwrap();

        println!("Seq: {},\nStart_fwd: {},\nEnd_fwd: {},\nHandles: {:#?}\n",
                 graph_kmer.seq, fwd_kmer.start, fwd_kmer.end, graph_kmer.handles.get(0).unwrap());
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

pub fn print_seq_to_file(seq : &str, path : &PathBuf) {
    let mut file = File::create(path).unwrap_or_else(|_| panic!("Error creating file {:?}", path));
    file.write_all(">ref\n".as_bytes());
    file.write_all(seq.as_bytes());
}

pub fn print_kmers_to_file(kmers : &Vec<Kmer>, path : &PathBuf) {
    let mut file = File::create(path).unwrap_or_else(|_| panic!("Error creating file {:?}", path));
    for i in 0..kmers.len() {
        let curr_kmer = kmers.get(i).unwrap();
        let index_str = format!(">kmer-{}\n", i);
        let kmers_str = format!("{}\n", curr_kmer.seq);
        file.write_all(index_str.as_bytes());
        file.write_all(kmers_str.as_bytes());
    }
}

pub fn print_kmers_to_file_split(kmers_on_fwd : &Vec<Kmer>, kmers_on_rev : &Vec<Kmer>, path : &PathBuf) {
    let mut file = File::create(path).unwrap_or_else(|_| panic!("Error creating file {:?}", path));
    for i in 0..kmers_on_fwd.len() {
        let curr_kmer = kmers_on_fwd.get(i).unwrap();
        let index_str = format!(">kmer-fwd-{}\n", i);
        let kmers_str = format!("{}\n", curr_kmer.seq);
        file.write_all(index_str.as_bytes());
        file.write_all(kmers_str.as_bytes());
    }
    for i in 0..kmers_on_rev.len() {
        let curr_kmer = kmers_on_rev.get(i).unwrap();
        let index_str = format!(">kmer-rev-{}\n", i);
        let kmers_str = format!("{}\n", curr_kmer.seq);
        file.write_all(index_str.as_bytes());
        file.write_all(kmers_str.as_bytes());
    }
}

pub fn verify_kmers(seq : &str, kmers : &Vec<Kmer>, path : &PathBuf) {
    let mut file = File::create(path).unwrap_or_else(|_| panic!("Error creating file {:?}", path));

    // This string will be used as reference to write the kmers in a way that makes sense
    let mut null_string = (0..seq.len()).map(|_| "-").collect::<String>();

    // Print reference to file
    let reference = format!("ref:       {}\n",seq);
    file.write_all(reference.as_bytes());

    // Print kmers to file
    let mut last_seq : usize = 0;
    let mut curr_kmer = kmers.get(last_seq).unwrap();
    for i in 0..seq.len()-3 {
        if seq[i..i+3] == curr_kmer.seq {
            let mut curr_string = null_string.clone();
            curr_string.replace_range((i..i+3), curr_kmer.seq.as_str());

            let kmer_string = format!("kmer:      {}\n",curr_string);
            file.write_all(kmer_string.as_bytes());

            last_seq += 1;
            if last_seq < kmers.len() {
                curr_kmer = kmers.get(last_seq).unwrap();
            } else {
                break;
            }
        }
    }
}

pub fn verify_kmers_2(seq : &str, kmers : &Vec<Kmer>, kmer_pos : &Vec<KmerPos>, path : &PathBuf) {
    let mut file = File::create(path).unwrap_or_else(|_| panic!("Error creating file {:?}", path));

    // This string will be used as reference to write the kmers in a way that makes sense
    let mut null_string = (0..seq.len()).map(|_| "-").collect::<String>();

    // Print reference to file
    let reference = format!("ref:       {}\n",seq);
    file.write_all(reference.as_bytes());

    for i in 0..kmers.len() {
        let curr_kmer = kmers.get(i).unwrap();
        let curr_pos = kmer_pos.get(i).unwrap();
        let mut curr_string = null_string.clone();

        curr_string.replace_range((curr_pos.start as usize..curr_pos.end as usize), curr_kmer.seq.as_str());
        let kmer_string = format!("kmer:      {}\n",curr_string);
        file.write_all(kmer_string.as_bytes());
    }

}