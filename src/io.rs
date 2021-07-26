use std::io::{Result, BufReader};
use std::io::prelude::*;
use std::fs::File;
use bstr::io::BufReadExt;
use substring::Substring;
use std::path::Path;
use std::ffi::OsStr;
use itertools::Itertools;
use rayon::prelude::*;
use std::io;

#[derive(PartialEq)]
enum InputFileTypes {
    Fasta,
    Fastq
}

pub struct InputSequence {
    pub name : String,
    pub seq : String
}

impl InputSequence {
    pub fn split_into_kmers(&self, kmer_size : usize) -> Vec<String> {
        let mut seq_kmers : Vec<String> = Vec::new();

        let query_string = String::from(self.seq.clone());
        let mut query_kmers : Vec<String> = Vec::new();
        for i in 0..(self.seq.len() - kmer_size)  {
            query_kmers.push(query_string.substring(i,i+kmer_size).into())
        }
        println!("Query kmers: {:#?}", seq_kmers);

        seq_kmers
    }
}

/// Parse a fasta/fastq file and returns the list of sequences from the given file
pub fn read_seqs_from_file(filename : &str) -> Result<Vec<InputSequence>> {
    let mut seqs : Vec<InputSequence> = Vec::new();

    let mut file = File::open(filename)?;
    let mut lines = BufReader::new(file).lines().peekable();

    // Check file type
    let file_type : InputFileTypes;
    let file_extension = Path::new(filename).extension().and_then(OsStr::to_str);
    match file_extension {
        Some("fasta") | Some("fa") => file_type = InputFileTypes::Fasta,
        Some("fastq") | Some("fq") => file_type = InputFileTypes::Fastq,
        _ => panic!("Unrecognized file type")
    }

    // Then parse the file itself
    if file_type == InputFileTypes::Fasta {
        while let (Some(Ok(name_long)), Some(Ok(seq))) = (lines.next(), lines.next()) {
                let name: String = String::from(
                    name_long.substring(1, name_long.len())
                );
                seqs.push(InputSequence { name, seq });
        }
    } else if file_type == InputFileTypes::Fastq {
        while let (Some(Ok(name)), Some(Ok(seq)), Some(Ok(_)), Some(Ok(_))) =
        (lines.next(), lines.next(), lines.next(), lines.next()) {
                seqs.push(InputSequence { name, seq });
        }
    }

    Ok(seqs)
}

#[cfg(test)]
mod test {
    use crate::io::read_seqs_from_file;

    #[test]
    fn test_read_fasta() {
        let test_seqs = read_seqs_from_file("./test/test.fa").unwrap();
        assert_eq!(test_seqs.len(), 10);
    }

    #[test]
    fn test_read_fastq() {
        let test_seqs = read_seqs_from_file("./test/test.fq").unwrap();
        assert_eq!(test_seqs.len(), 1);
    }

    /*
    #[test]
    fn test_wrong_extension() {
        let test_seqs = read_seqs_from_file("./test/test.gfa").unwrap();
        println!("Test: {:#?}", test_seqs);
    }
     */
}
