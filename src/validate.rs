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
        let nodes_ids = parse_nodes_from_path_matching(&alignment.path_matching.clone().unwrap());

        ValidationRecord {
            read_name: alignment.clone().query_name.unwrap(),
            CIGAR: alignment.clone().notes.unwrap(),
            read_seq: read.seq.clone().to_string(),
            nodes_id: nodes_ids.clone(),
            nodes_seq: get_nodes_sequences(graph, &nodes_ids),
        }
    }

    pub fn to_string(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{:#?}\n",
            self.read_name, self.CIGAR, self.read_seq, self.read_name, self.nodes_seq
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
            graph
                .sequence(Handle::from_integer(*id))
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
        //TODO: fix read to use
        .map(|a| ValidationRecord::from_graph_and_alignment(graph, a, reads.get(0).unwrap()))
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
}
