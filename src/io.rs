use std::ffi::OsStr;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, Result};
use std::path::Path;

use substring::Substring;

#[derive(PartialEq)]
enum InputFileTypes {
    Fasta,
    Fastq,
}

#[derive(Debug, Clone)]
pub struct QuerySequence {
    pub name: String,
    pub seq: String,
}

impl QuerySequence {
    pub fn new() -> Self {
        QuerySequence {
            name: "".to_string(),
            seq: "".to_string(),
        }
    }

    // TODO: maybe an Option return type would be better
    pub fn split_into_kmers(&self, kmer_size: usize) -> Vec<String> {
        let mut seq_kmers: Vec<String> = Vec::new();

        let query_string = self.seq.clone();

        // Check if it's possible to obtain kmers (otherwise return empty vec)
        if kmer_size <= query_string.len() {
            //let end = self.seq.len() - kmer_size;
            for i in 0..(self.seq.len() - kmer_size + 1) {
                seq_kmers.push(String::from(query_string.substring(i, i + kmer_size)));
            }
            //println!("Query kmers: {:#?}", seq_kmers);
        }

        seq_kmers
    }

    pub fn from_string(seq: &str) -> Self {
        QuerySequence {
            name: String::from(""),
            seq: seq.to_string().clone(),
        }
    }
}

/// Parse a fasta/fastq file and returns the list of sequences from the given file
pub fn read_seqs_from_file(filename: &str) -> Result<Vec<QuerySequence>> {
    let mut seqs: Vec<QuerySequence> = Vec::new();

    let file = File::open(filename)?;
    let mut lines = BufReader::new(file).lines().peekable();

    // Check file type
    let file_type: InputFileTypes;
    let file_extension = Path::new(filename).extension().and_then(OsStr::to_str);
    match file_extension {
        Some("fasta") | Some("fa") => file_type = InputFileTypes::Fasta,
        Some("fastq") | Some("fq") => file_type = InputFileTypes::Fastq,
        _ => panic!("Unrecognized file type"),
    }

    // Then parse the file itself
    if file_type == InputFileTypes::Fasta {
        /*
        while let (Some(Ok(name_long)), Some(Ok(seq))) = (lines.next(), lines.next()) {
            let name: String = String::from(name_long.substring(1, name_long.len()));
            seqs.push(QuerySequence { name, seq });
        }
         */
        let mut last_fasta_name = String::new();
        let mut fasta_same_name_count: i32 = 0;
        let mut last_fasta_seq = String::new();
        while let Some(Ok(line)) = lines.next() {
            if line.starts_with('>') {
                last_fasta_name = String::from(line.substring(1, line.len()));
                fasta_same_name_count = 0;
            } else {
                if line != "" {
                    last_fasta_seq = line;
                    match fasta_same_name_count == 0 {
                        true => seqs.push(QuerySequence {
                            name: last_fasta_name.clone(),
                            seq: last_fasta_seq,
                        }),
                        false => seqs.push(QuerySequence {
                            name: last_fasta_name.clone()
                                + fasta_same_name_count.to_string().as_str(),
                            seq: last_fasta_seq,
                        }),
                    }
                    fasta_same_name_count += 1;
                }
            }
        }
    } else if file_type == InputFileTypes::Fastq {
        while let (Some(Ok(name)), Some(Ok(seq)), Some(Ok(_)), Some(Ok(_))) =
            (lines.next(), lines.next(), lines.next(), lines.next())
        {
            seqs.push(QuerySequence { name, seq });
        }
    }

    /* <- this also reads fasta without headers, can be ignored
    if file_type == InputFileTypes::Fasta {
        while let Some(Ok(line)) = lines.next() {
            // TODO: keep or not? idk
            if line.starts_with(">") {
                let name: String = String::from(line.substring(1, line.len()));
                if let Some(Ok(seq)) = lines.next() {
                    seqs.push(QuerySequence { name, seq });
                }
            } else {
                // Just to be extra safe...
                if line != "" {
                    seqs.push(QuerySequence {
                        name: "".to_string(),
                        seq: line,
                    })
                }
            }
        }
    } else if file_type == InputFileTypes::Fastq {
        while let (Some(Ok(name)), Some(Ok(seq)), Some(Ok(_)), Some(Ok(_))) =
            (lines.next(), lines.next(), lines.next(), lines.next())
        {
            seqs.push(QuerySequence { name, seq });
        }
    }
     */

    Ok(seqs)
}

#[cfg(test)]
mod test {
    use crate::io::{read_seqs_from_file, QuerySequence};

    #[test]
    fn test_read_fasta_single_read() {
        let test_seqs = read_seqs_from_file("./test/single-read-test.fa").unwrap();
        println!("Test seqs: {:#?}", test_seqs);
        assert_eq!(test_seqs.len(), 1);
        assert_eq!(test_seqs.get(0).unwrap().name, "seq0".to_string());
        assert_eq!(
            test_seqs.get(0).unwrap().seq,
            "AAAAACGTTAAATTTGGCATCGTAGCAAAAA"
        );
    }

    #[test]
    fn test_read_fasta_headers() {
        let test_seqs = read_seqs_from_file("./test/multiple-read-test.fa").unwrap();
        assert_eq!(test_seqs.len(), 2);
        assert_eq!(test_seqs.get(0).unwrap().name, "seq0".to_string());
        assert_eq!(
            test_seqs.get(0).unwrap().seq,
            "AAAAACGTTAAATTTGGCATCGTAGCAAAAA"
        );
        assert_eq!(test_seqs.get(1).unwrap().name, "seq1".to_string());
        assert_eq!(test_seqs.get(1).unwrap().seq, "TTTCGTTAAATTTGGCATCGTAGCTTT");
    }

    /*
    #[test]
    fn test_read_fasta_no_header() {
        let test_seqs = read_seqs_from_file("./test/test-no-headers.fa").unwrap();
        assert_eq!(test_seqs.len(), 2);
        assert_eq!(test_seqs.get(0).unwrap().name, "".to_string());
        assert_eq!(
            test_seqs.get(0).unwrap().seq,
            "AAAAACGTTAAATTTGGCATCGTAGCAAAAA"
        );
        assert_eq!(test_seqs.get(1).unwrap().name, "".to_string());
        assert_eq!(test_seqs.get(1).unwrap().seq, "TTTCGTTAAATTTGGCATCGTAGCTTT");
    }
     */

    #[test]
    fn test_read_fastq() {
        let test_seqs = read_seqs_from_file("./test/test.fq").unwrap();
        assert_eq!(test_seqs.len(), 1);
    }

    #[test]
    fn test_split_ok() {
        let test_seq = QuerySequence::from_string(&String::from("AAACTG"));
        let seq_kmers = test_seq.split_into_kmers(3);
        assert_eq!(seq_kmers.len(), 4);
        assert_eq!(*seq_kmers.get(0).unwrap(), String::from("AAA"));
        assert_eq!(*seq_kmers.get(1).unwrap(), String::from("AAC"));
        assert_eq!(*seq_kmers.get(2).unwrap(), String::from("ACT"));
        assert_eq!(*seq_kmers.get(3).unwrap(), String::from("CTG"));
    }

    #[test]
    fn test_split_greater() {
        let test_seq = QuerySequence::from_string(&String::from("AAA"));
        let seq_kmers = test_seq.split_into_kmers(4);
        assert_eq!(seq_kmers.len(), 0);
    }

    #[test]
    fn test_split_lesser() {
        let test_seq = QuerySequence::from_string(&String::from("AA"));
        let seq_kmers = test_seq.split_into_kmers(3);
        assert_eq!(seq_kmers.len(), 0);
    }
}
