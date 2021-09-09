use crate::chain::Chain;
use crate::index::Index;
use crate::io::QuerySequence;
use crate::kmer::{SeqOrient, SeqPos};
use ab_poa::abpoa_wrapper::{AbpoaAligner, AbpoaAlignmentResult};
use bstr::ByteVec;
use core::cmp;
use handlegraph::handle::Handle;
use handlegraph::handlegraph::HandleGraph;
use handlegraph::hashgraph::HashGraph;
use itertools::Itertools;
use rayon::prelude::*;
use std::ops::Range;
use substring::Substring;

/*
fn find_range_single_chain(index: &Index, chain: &Chain) -> Range<u64> {
    let min_node = chain
        .anchors
        .par_iter()
        .map(|a| index.node_id_from_seqpos(&a.target_begin))
        .min()
        .unwrap();
    let max_node = chain
        .anchors
        .par_iter()
        .map(|a| index.node_id_from_seqpos(&a.target_end))
        .max()
        .unwrap();
    min_node..max_node
}
pub fn extract_query_subsequence(
    index: &Index,
    chains: &Vec<Chain>,
    _query: &QuerySequence,
) -> Vec<String> {
    let subseqs: Vec<String> = chains
        .par_iter()
        .map(|c| index.seq_from_start_end_seqpos(&c.target_begin, &c.target_end))
        .collect();
    subseqs
}
 */

pub(crate) fn obtain_base_level_alignment(index: &Index, chain: &Chain) -> GAFAlignment {
    // Find the range of node ids involved in the alignment
    let po_range = find_range_chain(index, chain);

    // Find nodes and edges
    let (nodes, edges) = find_nodes_edges_for_abpoa(&index, &po_range);
    // TODO: possibly avoid this
    let nodes_str: Vec<&str> = nodes.par_iter().map(|x| &x as &str).collect();

    println!("Seqs: {:?}", nodes_str);
    println!("Edges: {:?}", edges);

    // Find subquery implied by the chain
    let subquery_range = chain.find_query_start_end();
    let subquery = chain
        .query
        .seq
        .to_string()
        .substring(subquery_range.start as usize, subquery_range.end as usize)
        .to_string();
    println!("Subquery is: {:#?}", subquery);

    // Align with abpoa
    let mut result = AbpoaAlignmentResult::new();
    unsafe {
        result = align_with_poa(&nodes_str, &edges, subquery.as_str());
    }
    println!("Result is: {:#?}", result);
    let alignment: GAFAlignment = generate_alignment(
        chain,
        &result,
        &po_range,
        &subquery_range,
        chain.query.seq.len(),
    );

    alignment
}

/*
fn find_range_single_chain_anchor(index: &Index, chain: &Chain) -> Range<u64> {
    let min_handle = chain
        .anchors
        .par_iter()
        .map(|a| index.handle_from_seqpos(&a.target_begin))
        .min()
        .unwrap();
    let max_handle = chain
        .anchors
        .par_iter()
        .map(|a| index.handle_from_seqpos(&a.target_end))
        .max()
        .unwrap();
    u64::from(min_handle)..u64::from(max_handle)
}
 */

enum RangeOrient {
    Forward,
    Reverse,
    Both,
}

struct OrientedGraphRange {
    pub orient: RangeOrient,
    pub handles: Vec<Handle>,
}
impl OrientedGraphRange {
    pub fn get_first_handle(&self) -> &Handle {
        self.handles.get(0).unwrap()
    }
    pub fn get_last_handle(&self) -> &Handle {
        self.handles.last().unwrap()
    }
}

fn find_range_chain(index: &Index, chain: &Chain) -> OrientedGraphRange {
    let min_handle: Handle = chain
        .anchors
        .par_iter()
        .map(|a| index.handle_from_seqpos(&a.target_begin))
        .min()
        .unwrap();
    let max_handle: Handle = chain
        .anchors
        .par_iter()
        .map(|a| index.handle_from_seqpos(&a.target_end))
        .max()
        .unwrap();

    let mut po_range_handles: Vec<Handle> = Vec::new();
    let mut orient: RangeOrient = RangeOrient::Both;

    /*
    println!(
        "Min handle: {:#?}, Max handle: {:#?}",
        min_handle, max_handle
    );
    println!("Range as handle is: {:#?}", min_handle..max_handle);
    println!(
        "Range as nodes is: {:#?}",
        u64::from(min_handle)..u64::from(max_handle)
    );
     */

    if !min_handle.is_reverse() && !max_handle.is_reverse() {
        po_range_handles = (u64::from(min_handle)..u64::from(max_handle))
            .into_iter()
            .filter(|x| x % 2 == 0)
            .map(|x| Handle::from_integer(x * 2))
            .collect();
        orient = RangeOrient::Forward;
    } else if min_handle.is_reverse() && max_handle.is_reverse() {
        po_range_handles = (u64::from(min_handle)..u64::from(max_handle))
            .into_iter()
            .filter(|x| x % 2 != 0)
            .map(|x| Handle::from_integer(x * 2 + 1))
            .collect();
        orient = RangeOrient::Reverse;
    } else {
        let po_range_handles_fwd: Vec<Handle> = (u64::from(min_handle)..u64::from(max_handle))
            .into_iter()
            .map(|x| Handle::from_integer(x * 2))
            .collect();
        let po_range_handles_rev: Vec<Handle> = (u64::from(min_handle)..u64::from(max_handle))
            .into_iter()
            .map(|x| Handle::from_integer(x * 2 + 1))
            .collect();
        po_range_handles = [po_range_handles_fwd, po_range_handles_rev].concat();
        po_range_handles.sort();

        orient = RangeOrient::Both
    }

    OrientedGraphRange {
        orient,
        handles: po_range_handles,
    }
}

/*
/// Finds the nodes and edges involved in the alignment
fn find_nodes_edges_2_anchors(
    index: &Index,
    po_range: &Range<u64>,
) -> (Vec<String>, Vec<(u64, u64)>) {
    let seqs: Vec<String> = po_range
        .clone()
        .into_iter()
        .map(|h| index.seq_from_handle(&Handle::from_integer(h)))
        .collect();
    println!("Seqs: {:#?}", seqs);

    let edges: Vec<(u64, u64)> = po_range
        .clone()
        .into_iter()
        .flat_map(|handle| {
            index
                .outgoing_edges_from_handle(&Handle::from_integer(handle))
                .par_iter()
                .map(|x| x.as_integer())
                .filter(|target_handle| po_range.contains(target_handle))
                .map(move |target_handle| (handle - po_range.start, target_handle - po_range.start))
                .collect::<Vec<(u64, u64)>>()
        })
        .collect();

    (seqs, edges)
}
 */

/// Finds the nodes and edges involved in the alignment
fn find_nodes_edges_for_abpoa(
    index: &Index,
    po_range: &OrientedGraphRange,
) -> (Vec<String>, Vec<(usize, usize)>) {
    let range_handles = po_range.handles.clone();

    // Get the seqs for the handles in the po_range
    // This is necessary because in abpoa the nodes are specified by a Vec of sequences
    let seqs: Vec<String> = range_handles
        .par_iter()
        .map(|h| index.seq_from_handle(h))
        .collect();

    // Get the edges between the handles in the po_range
    // In abpoa the nodes are given as 0-based tuples (starting_node, ending_node).
    // For example in order to create this graph: "ACG" -> "AAA"
    // the following instuction will be used:
    // aligner.add_nodes_edges(&vec!["ACG", "AAA"], &vec![(0, 1)]);
    let edges: Vec<(usize, usize)> = range_handles
        .clone()
        .into_iter()
        // For each handle in the po range
        .flat_map(|handle| {
            index
                // Find all the outgoing edges (this actually returns the handles the edges point to)
                .outgoing_edges_from_handle(&handle)
                .into_iter()
                // Only keep the ones in the po range
                .filter(|target_handle| range_handles.contains(target_handle))
                // And find their 0-based position in the po range
                .map(|target_handle| {
                    let starting_pos = range_handles.par_iter().position(|x| *x == handle).unwrap();
                    let ending_pos = range_handles
                        .par_iter()
                        .position(|x| *x == target_handle)
                        .unwrap();
                    (starting_pos, ending_pos)
                }) //(handle.as_integer() - po_range.par_iter().min().unwrap().as_integer(), target_handle.as_integer() - po_range.par_iter().min().unwrap().as_integer()))
                .collect::<Vec<(usize, usize)>>()
        })
        .collect();

    (seqs, edges)
}

/*
/// Finds the nodes and edges involved in the alignment
fn find_nodes_edges_2_anchors_vec_incoming(
    index: &Index,
    po_range: &Vec<Handle>,
) -> (Vec<String>, Vec<(usize, usize)>) {
    // Get the seqs for the handles in the po_range
    // This is necessary because in abpoa the nodes are specified by a Vec of sequences
    let seqs: Vec<String> = po_range.par_iter().map(|h| index.seq_from_handle(h)).collect();

    // Get the edges between the handles in the po_range
    // In abpoa the nodes are given as 0-based tuples (starting_node, ending_node).
    // For example in order to create this graph: "ACG" -> "AAA"
    // the following instuction will be used:
    // aligner.add_nodes_edges(&vec!["ACG", "AAA"], &vec![(0, 1)]);
    let edges: Vec<(usize, usize)> = po_range
        .clone()
        .into_iter()
        // For each handle in the po range
        .flat_map(|handle| {
            index
                // Find all the outgoing edges (this actually returns the handles the edges point to)
                .incoming_edges_from_handle(&handle)
                .into_iter()
                // Only keep the ones in the po range
                .filter(|target_handle| po_range.contains(target_handle))
                // And find their 0-based position in the po range
                .map(move |target_handle| {
                    let ending_pos = po_range.par_iter().position(|x| *x == handle).unwrap();
                    let starting_pos = po_range.par_iter().position(|x| *x == target_handle).unwrap();
                    (starting_pos, ending_pos)
                }) //(handle.as_integer() - po_range.par_iter().min().unwrap().as_integer(), target_handle.as_integer() - po_range.par_iter().min().unwrap().as_integer()))
                .collect::<Vec<(usize, usize)>>()
        })
        .collect();

    (seqs, edges)
}
 */

/*
    GAF format
    https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf
    Col     Type    Description
    1       string  Query sequence name
    2       int     Query sequence length
    3       int     Query start (0-based; closed)
    4       int     Query end (0-based; open)
    5       char    Strand relative to the path: "+" or "-"
    6       string  Path matching /([><][^\s><]+(:\d+-\d+)?)+|([^\s><]+)/
    7       int     Path length
    8       int     Start position on the path (0-based)
    9       int     End position on the path (0-based)
    10      int     Number of residue matches
    11      int     Alignment block length
    12      int     Mapping quality (0-255; 255 for missing)
*/
#[derive(Debug)]
pub(crate) struct GAFAlignment {
    query_name: String,
    query_length: u64,
    query_start: u64,
    query_end: u64,
    strand: char,
    path_matching: String,
    path_length: u64,
    path_start: u64,
    path_end: u64,
    residue: u64,
    alignment_block_length: u64,
    mapping_quality: u64,
    notes: String,
}
impl GAFAlignment {
    pub(crate) fn from_chain(chain: &Chain) -> Self {
        GAFAlignment {
            query_name: chain.query.name.clone(),
            query_length: chain.query.seq.len() as u64,
            query_start: chain.anchors.front().unwrap().query_begin,
            query_end: chain.anchors.back().unwrap().query_end,
            strand: '+',
            path_matching: " ".to_string(),
            path_length: 0,
            path_start: 0,
            path_end: 0,
            residue: 0,
            alignment_block_length: 0,
            mapping_quality: cmp::min(chain.mapping_quality as u64, 254_u64),
            notes: "ta:Z:chain".to_string(),
        }
    }

    pub fn to_string(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            self.query_name,
            self.query_length,
            self.query_start,
            self.query_end,
            self.strand,
            self.path_matching,
            self.path_length,
            self.path_start,
            self.path_end,
            self.residue,
            self.alignment_block_length,
            self.mapping_quality,
            self.notes
        )
    }
}

/*
fn find_nodes_edges(graph: &HashGraph, po_range: &Range<u64>) -> (Vec<String>, Vec<(u64, u64)>) {
    let mut handles: Vec<Handle> = graph.handles_iter().collect();
    handles.par_sort();
    let mut node_handles: Vec<Handle> = Vec::new();
    for val in po_range.start..po_range.end {
        node_handles.push(*handles.get((val - 1) as usize).unwrap());
    }
    let mut node_seqs: Vec<String> = Vec::new();
    for handle in &node_handles {
        node_seqs.push(graph.sequence(*handle).into_string_lossy());
    }
    let mut edges: Vec<(u64, u64)> = node_handles
        .par_iter()
        .combinations(2)
        .filter_map(|h| {
            let left = **h.get(0).unwrap();
            let right = **h.get(1).unwrap();
            if graph.has_edge(left, right) {
                Some((u64::from(left), u64::from(left)))
            } else {
                None
            }
        })
        .collect();
    // Since in abpoa nodes start from 0, the edges have to be "normalized".
    // In practice this means subtracting from each each edge the start of the po_range,
    // e.g. po_range: [4,7), edges: [(4,5), (5,6), (6,7)] -> norm_edges [(0,1), (1,2), (2, 3)]
    let norm_edges: Vec<(u64, u64)> = edges
        .into_par_iter()
        .map(|x| {
            let start = x.0;
            let end = x.1;
            (start - po_range.start, end - po_range.start)
        })
        .collect();
    (node_seqs, norm_edges)
}
*/

/// Use rs-abpoa to align the sequences to the graph
unsafe fn align_with_poa(
    nodes: &Vec<&str>,
    edges: &Vec<(usize, usize)>,
    query: &str,
) -> AbpoaAlignmentResult {
    let mut aligner = AbpoaAligner::new_with_example_params();
    aligner.add_nodes_edges(nodes, edges);
    let res = aligner.align_sequence(query);
    res
}
// TODO: complete this
// The final alignment should be made up by:
// - perfect matches (= anchors) where applicable
// - poa result (cigar) otherwise
fn generate_alignment(
    chain: &Chain,
    result: &AbpoaAlignmentResult,
    po_range: &OrientedGraphRange,
    subquery_range: &Range<u64>,
    og_query_length: usize,
) -> GAFAlignment {
    // Get the nodes involved in the alignment from abPOA
    let mut alignment_path = result.graph_nodes.clone();
    // Since abpoa deals with 1-base nodes, the abstraction simply converts
    // the abpoa_node_id to their position in the abstraction
    // i.e suppose the abpoa_nodes are 2,3,4,5; the abstraction nodes are [2] and [3,4,5];
    // result.graph_nodes will contain [0, 1, 1, 1], I want: [0,1]
    alignment_path.dedup();

    // Since alignment_path contains ids that are coherent with the nodes/edges
    // used during its creation, the nodes will be something like [0, 1, ... po_range.max].
    // However, this graph is a subgraph of our original graph, so the ids must be shifted
    // accordingly.
    let og_graph_alignment_path: Vec<Handle> = alignment_path
        .par_iter()
        .map(|x| *po_range.handles.get(*x).unwrap())
        .collect();

    let alignment_path_string: Vec<String> = og_graph_alignment_path
        .par_iter()
        .map(|handle| match handle.is_reverse() {
            false => ">".to_string() + handle.as_integer().to_string().as_str(),
            true => "<".to_string() + handle.as_integer().to_string().as_str(),
        })
        .collect();

    // Find path start and end
    let og_path_start = og_graph_alignment_path.first().unwrap().as_integer();
    let og_path_end = og_graph_alignment_path.last().unwrap().as_integer();

    GAFAlignment {
        query_name: chain.query.name.clone(),
        query_length: og_query_length as u64,
        query_start: subquery_range.start,
        query_end: subquery_range.end,
        strand: '+',
        path_matching: alignment_path_string.iter().join(""),
        path_length: og_path_end - og_path_start,
        path_start: og_path_start,
        path_end: og_path_end,
        residue: 0,
        alignment_block_length: result.n_aligned_bases as u64,
        mapping_quality: 255, //result.best_score as u64,
        notes: ("cg:Z:".to_string()
            + result.cigar.clone().as_str()
            + ":as:i:"
            + result.best_score.to_string().as_str()),
    }
}

#[cfg(test)]
mod test {
    use crate::chain::Chain;
    use crate::index::Index;
    use gfa::gfa::GFA;
    use gfa::parser::GFAParser;
    use handlegraph::hashgraph::HashGraph;

    /*
    fn graph_from_file() -> HashGraph {
        let parser = GFAParser::new();
        let gfa: GFA<usize, ()> = parser
            .parse_file(&PathBuf::from("./test/test.gfa"))
            .unwrap();
        let graph = HashGraph::from_gfa(&gfa);
        graph
    }

    #[test]
    fn test_po_range() {
        let graph = graph_from_file();
        let index = Index::build(&graph, 11, 100, 100, 7.0, None);
    }

     */
}
