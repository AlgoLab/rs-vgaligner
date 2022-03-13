use crate::align::GAFAlignment;
use crate::index::Index;
use crate::io::QuerySequence;
use bstr::ByteVec;
use handlegraph::handle::Handle;
use handlegraph::handlegraph::HandleGraph;
use handlegraph::hashgraph::HashGraph;
use itertools::Itertools;
use regex::Regex;
use std::fs::File;
use std::io::Write;
use log::info;

pub struct ValidationRecord {
    pub read_name: String,
    pub CIGAR: String,
    pub read_seq: String,
    pub nodes_id: Vec<u64>,
    pub nodes_seq: Vec<String>,
}

impl ValidationRecord {
    pub fn new() -> Self {
        ValidationRecord {
            read_name: "".to_string(),
            CIGAR: "".to_string(),
            read_seq: "".to_string(),
            nodes_id: vec![],
            nodes_seq: vec![],
        }
    }

    pub fn from_graph_and_alignment(
        graph: &HashGraph,
        alignment: &GAFAlignment,
        read: &QuerySequence,
    ) -> Self {
        info!("Validating alignment for {}", alignment.clone().query_name.unwrap());

        match alignment.clone().path_matching {
            Some(path) => {
                //info!("Path matching {}", path);
                let nodes_ids = parse_nodes_from_path_matching(&path);

                ValidationRecord {
                    read_name: alignment.clone().query_name.unwrap(),
                    CIGAR: alignment.clone().notes.unwrap().split(',').last().unwrap().to_string(),
                    read_seq: read.seq.clone().to_string(),
                    nodes_id: nodes_ids.clone(),
                    nodes_seq: get_nodes_sequences(graph, &nodes_ids),
                }
            },
            _ => {
                info!("Read {} could not be aligned", alignment.clone().query_name.unwrap());

                ValidationRecord {
                    read_name: alignment.clone().query_name.unwrap(),
                    CIGAR: "NOT ALIGNED".to_string(),
                    read_seq: read.seq.clone().to_string(),
                    nodes_id: vec![],
                    nodes_seq: vec![]
                }
            }
        }

    }

    pub fn to_string(&self) -> String {
        /*
        format!(
            "{}\t{}\t{}\t{:#?}\t{:#?}\n",
            self.read_name, self.CIGAR, self.read_seq, self.nodes_id, self.nodes_seq
        )
         */

        format!(
            "{}\n{}\n{}\n{:?}\n{:?}\n\n",
            self.read_name, self.CIGAR, self.read_seq, self.nodes_id, self.nodes_seq
        )
    }
}

pub fn parse_nodes_from_path_matching(path_matching: &str) -> Vec<u64> {
    let re = Regex::new(r"(>|<)([0-9]+)").unwrap();

    re.captures_iter(path_matching)
        .map(|node| node[2].parse::<u64>().unwrap())
        .collect()
}

pub fn get_nodes_sequences(graph: &HashGraph, nodes_ids: &Vec<u64>) -> Vec<String> {
    nodes_ids
        .iter()
        .map(|id| {
            //info!("id is: {}", id);
            //TODO: validation should take strand into account,
            //not necessary but would probably be better
            graph
                .sequence(Handle::pack(*id, false))
                .into_string_lossy()
        })
        .collect()
}

pub fn create_validation_records(
    graph: &HashGraph,
    alignments: &Vec<GAFAlignment>,
    reads: &Vec<QuerySequence>,
) -> Vec<ValidationRecord> {
    let records: Vec<ValidationRecord> = alignments
        .iter()
        .map(|a| {
            let aligned_read = reads.iter().find(|x| x.name == a.clone().query_name.unwrap()).unwrap();
            ValidationRecord::from_graph_and_alignment(graph, a, aligned_read)
        })
        .collect();
    records
}

pub fn write_validation_to_file(
    validation_records: &Vec<ValidationRecord>,
    file_name: String,
) -> std::io::Result<()> {
    let val_strings: Vec<String> = validation_records
        .iter()
        .map(|val| val.to_string())
        .collect();
    let mut file =
        File::create(&file_name).unwrap_or_else(|_| panic!("Couldn't create file {}", &file_name));
    file.write_all(&val_strings.join("").as_bytes())
        .unwrap_or_else(|_| panic!("Couldn't write to file {}", &file_name));
    Ok(())
}

#[cfg(test)]
mod test {
    use crate::validate::parse_nodes_from_path_matching;

    #[test]
    fn test_simple_parsing() {
        assert_eq!(parse_nodes_from_path_matching(">1<2>3"), vec![1, 2, 3])
    }

    #[test]
    fn test_double_digit_parsing() {
        assert_eq!(parse_nodes_from_path_matching(">10<20"), vec![10, 20])
    }

    #[test]
    fn test_empty_parsing() { assert_eq!(parse_nodes_from_path_matching("*"),  Vec::<u64>::new()) }
}
