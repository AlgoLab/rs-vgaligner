use crate::dna::reverse_complement;
use crate::kmer::{
    generate_hash, generate_kmers, generate_kmers_linearly,
    generate_pos_on_ref, Kmer, KmerPos,
};
use crate::utils::{find_graph_seq_length, find_sequence_po, NodeRef};
use ahash::RandomState;
use boomphf::hashmap::BoomHashMap;
use bv::BitVec;
use handlegraph::handle::{Edge, Handle};
use handlegraph::hashgraph::HashGraph;
use crate::io::store_object_to_file;
use serde_with::serde_as;
use boomphf::hashmap::NoKeyBoomHashMap;
//use crate::serialization::KmerTable;

#[derive(Default)]
pub struct Index {
    // the kmer size that this graph was built on
    kmer_length: u64,
    // consider only kmers where to_key(kmer) % sampling_mod == 0
    sampling_length: u64,
    // total sequence length of the graph
    seq_length: u64,
    // forward sequence of the graph, stored here for fast access during alignment
    seq_fwd: String,
    // reverse complemented sequence of the graph, for fast access during alignment
    seq_rev: String,
    // mark node starts
    seq_bv: BitVec,
    // lets us map between our seq vector and handles (when our input graph is compacted!)
    //seq_by_rank: BitVec, //TODO: check this
    // edge count
    n_edges: u64,
    // what's the most-efficient graph topology we can store?
    edges: Vec<Edge>,
    // node count
    n_nodes: u64,
    // refer to ranges in edges
    node_ref: Vec<NodeRef>,
    // number of kmers in the index
    n_kmers: u64,
    // number of kmer positions in the index
    n_kmer_pos: u64,
    // our kmer hash table
    //bhpf :  TODO: hash
    // our kmer reference table (maps from bphf to index in kmer_pos_vec)
    //kmer_pos_ref: Vec<u64>,
    // our kmer positions
    kmer_pos_table: Vec<KmerPos>,
    // if we're loaded, helps during teardown
    loaded: bool,
}

impl Index {
    //See this pattern here: https://stackoverflow.com/a/41510505/5627359
    fn new() -> Self {
        Default::default()
    }

    pub fn build(
        graph: &HashGraph,
        kmer_length: u64,
        max_furcations: u64,
        max_degree: u64,
        sampling_rate: f32,
        out_prefix: &str,
    ) -> Self {
        // Get the number of nodes in the graph
        let number_nodes = graph.graph.len();

        // Get the length of the sequence encoded by the graph
        let seq_length = find_graph_seq_length(graph);

        // Mark node starts in forward
        let mut seq_bv: BitVec = BitVec::new_fill(false, seq_length + 1);

        // Store offsets in fwd and edge vector
        let mut node_ref: Vec<NodeRef> = Vec::new();

        let seq_fwd = find_sequence_po(graph, &mut seq_bv, &mut node_ref);
        let seq_rev = reverse_complement(&seq_fwd.as_str());

        // Store forward and reverse to file
        let seq_fwd_filename : String = out_prefix.to_owned() + ".sqf";
        store_object_to_file(seq_fwd.as_bytes(), &seq_fwd_filename);
        let seq_rev_filename : String = out_prefix.to_owned() + ".sqr";
        store_object_to_file(seq_rev.as_bytes(), &seq_rev_filename);

        // Store noderef to file
        let node_ref_filename : String = out_prefix.to_owned() + ".gyn";
        let encoded_node_ref = bincode::serialize(&node_ref).unwrap();
        store_object_to_file(&encoded_node_ref, &node_ref_filename);

        // Store seq_bv to file
        let bitvec_filename : String = out_prefix.to_owned() + ".sbv";
        let encoded_bit_vec = bincode::serialize(&seq_bv).unwrap();
        store_object_to_file(&encoded_bit_vec, &bitvec_filename);

        let kmers_on_graph: Vec<Kmer>;

        if graph.paths.is_empty() {
            // No info on paths available
            kmers_on_graph = generate_kmers(
                graph,
                kmer_length as u64,
                Some(max_furcations),
                Some(max_degree),
            );
        } else {
            // Optimize by using paths
            kmers_on_graph = generate_kmers_linearly(
                graph,
                kmer_length as u64,
                Some(max_furcations),
                Some(max_degree),
            );
        }

        // Store graph kmers to file
        let kmer_graph_filename : String = out_prefix.to_owned() + ".kgph";
        let encoded_kmer_graph = bincode::serialize(&kmers_on_graph).unwrap();
        store_object_to_file(&encoded_kmer_graph, &kmer_graph_filename);


        let kmers_positions_on_ref: Vec<KmerPos> =
            generate_pos_on_ref(&graph, &kmers_on_graph, &seq_length, &node_ref);

        // Store kmers_positions_on_ref to file
        let kmer_pos_filename : String = out_prefix.to_owned() + ".kpos";
        let encoded_kmer_pos = bincode::serialize(&kmers_positions_on_ref).unwrap();
        store_object_to_file(&encoded_kmer_pos, &kmer_pos_filename);

        // First obtain the kmers' hashes
        let hash_build = RandomState::with_seeds(0,0,0,0);
        let hashes = generate_hash(&kmers_on_graph, &hash_build);
        println!("{:#?}", hashes);

        // Then generate the table
        let kmers_table = NoKeyBoomHashMap::new_parallel(hashes.clone(), kmers_positions_on_ref.clone());
        println!("{:#?}", kmers_table);

        // Store the table to file
        let table_filename : String = out_prefix.to_owned() + ".bbx";
        let encoded_table = bincode::serialize(&kmers_table).unwrap();
        store_object_to_file(&encoded_table, &table_filename);

        Index {
            kmer_length,
            sampling_length: 0,
            seq_length,
            seq_fwd,
            seq_rev,
            seq_bv,
            //seq_by_rank: Default::default(),
            n_edges: 0,
            edges: vec![],
            n_nodes: 0,
            node_ref,
            n_kmers: kmers_positions_on_ref.len() as u64,
            n_kmer_pos: 0,
            //kmer_pos_ref: vec![],
            kmer_pos_table: kmers_positions_on_ref,
            loaded: true
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::kmer::{generate_kmers_linearly, generate_kmers_linearly_2};
    use handlegraph::hashgraph::Node;
    use handlegraph::mutablehandlegraph::MutableHandleGraph;
    use handlegraph::pathgraph::PathHandleGraph;
    use substring::Substring;

    /// This function creates a simple graph, used for debugging
    ///        | 2: CT \
    /// 1:  A            4: GCA
    ///        \ 3: GA |
    fn create_simple_graph() -> HashGraph {
        let mut graph: HashGraph = HashGraph::new();

        let h1 = graph.create_handle("A".as_bytes(), 1);
        let h2 = graph.create_handle("CT".as_bytes(), 2);
        let h3 = graph.create_handle("GA".as_bytes(), 3);
        let h4 = graph.create_handle("GCA".as_bytes(), 4);

        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h1, h3));
        graph.create_edge(&Edge(h2, h4));
        graph.create_edge(&Edge(h3, h4));

        let p1 = graph.create_path_handle("P1".as_bytes(), false);
        graph.append_step(&p1, h1);
        graph.append_step(&p1, h2);
        graph.append_step(&p1, h4);

        let p2 = graph.create_path_handle("P2".as_bytes(), false);
        graph.append_step(&p2, h1);
        graph.append_step(&p2, h3);
        graph.append_step(&p2, h4);

        graph
    }

    /// This function creates a simple graph, used for debugging
    ///        | 2: T \
    /// 1:  GAT            4: CA
    ///        \ 3: A |
    fn create_simple_graph_2() -> HashGraph {
        let mut graph: HashGraph = HashGraph::new();
        let h1 = graph.create_handle("GAT".as_bytes(), 1);
        let h2 = graph.create_handle("T".as_bytes(), 2);
        let h3 = graph.create_handle("A".as_bytes(), 3);
        let h4 = graph.create_handle("CA".as_bytes(), 4);

        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h1, h3));
        graph.create_edge(&Edge(h2, h4));
        graph.create_edge(&Edge(h3, h4));

        graph
    }

    #[test]
    fn test_simple_graph_2() {
        let graph = create_simple_graph_2();

        let total_length = find_graph_seq_length(&graph);
        let mut seq_bv: BitVec = BitVec::new_fill(false, total_length + 1);
        let mut node_ref: Vec<NodeRef> = Vec::new();

        let forward = find_sequence_po(&graph, &mut seq_bv, &mut node_ref);

        let kmers_graph = generate_kmers(&graph, 3, Some(100), Some(100));
        let kmers_ref = generate_pos_on_ref(&graph, &kmers_graph, &total_length, &node_ref);

        assert_eq!(kmers_graph.len(), kmers_ref.len());

        for i in 0..kmers_graph.len() {
            let graph_kmer = kmers_graph.get(i).unwrap();
            let ref_kmer = kmers_ref.get(i).unwrap();

            println!("{:#?}", graph_kmer);
            println!("{:#?}", ref_kmer);
            println!();
        }
    }

    #[test]
    fn test_assert_both_functions_find_same_kmers() {
        let graph = create_simple_graph();

        let total_length = find_graph_seq_length(&graph);
        let mut seq_bv: BitVec = BitVec::new_fill(false, total_length + 1);
        let mut node_ref: Vec<NodeRef> = Vec::new();

        let forward = find_sequence_po(&graph, &mut seq_bv, &mut node_ref);

        let kmers_graph_dozyg = generate_kmers(&graph, 3, Some(100), Some(100));
        let kmers_graph_rust_ver = generate_kmers_linearly(&graph, 3, Some(100), Some(100));

        assert_eq!(kmers_graph_dozyg.len(), kmers_graph_rust_ver.len());

        for kmer in &kmers_graph_rust_ver {
            assert!(kmers_graph_dozyg.contains(kmer));
        }
    }

    #[test]
    fn test_forward_creation() {
        let graph = create_simple_graph();

        let total_length = find_graph_seq_length(&graph);
        let mut seq_bv: BitVec = BitVec::new_fill(false, total_length + 1);
        let mut node_ref: Vec<NodeRef> = Vec::new();

        assert_eq!(total_length, 8);
        assert_eq!(
            "ACTGAGCA",
            find_sequence_po(&graph, &mut seq_bv, &mut node_ref)
        );

        use bv::*;
        // The bitvector marks node starts, so since nodes in the linearization are
        // A -> CT -> GA -> GCA the bitvector should be...
        assert_eq!(
            bit_vec![true, true, false, true, false, true, false, false, true],
            seq_bv
        );

        assert_eq!(node_ref.len(), 4);

        assert_eq!(
            *node_ref.get(0).unwrap(),
            NodeRef {
                seq_idx: 0,
                edge_idx: 0,
                edges_to_node: 0
            }
        );
        assert_eq!(
            *node_ref.get(1).unwrap(),
            NodeRef {
                seq_idx: 1,
                edge_idx: 2,
                edges_to_node: 1
            }
        );
        assert_eq!(
            *node_ref.get(2).unwrap(),
            NodeRef {
                seq_idx: 3,
                edge_idx: 3,
                edges_to_node: 1
            }
        );
        assert_eq!(
            *node_ref.get(3).unwrap(),
            NodeRef {
                seq_idx: 5,
                edge_idx: 4,
                edges_to_node: 2
            }
        );
    }

    #[test]
    fn test_kmers_graph_generation() {
        let graph = create_simple_graph();

        let kmers_on_graph = generate_kmers(&graph, 3, Some(100), Some(100));

        assert_eq!(kmers_on_graph.len(), 12);

        let kmers_on_graph_100 = generate_kmers(&graph, 6, Some(100), Some(100));
        assert_eq!(kmers_on_graph_100.len(), 4);

        // Check if it crashes with k=100...
        let kmers_on_graph_100 = generate_kmers(&graph, 100, Some(100), Some(100));
        assert_eq!(kmers_on_graph_100.len(), 0);
        // ...it doesn't!
    }

    #[test]
    fn test_simple_path() {
        let mut graph = HashGraph::new();

        // ACG -> TTT -> CA
        let h1 = graph.create_handle("ACG".as_bytes(), 1);
        let h2 = graph.create_handle("TTT".as_bytes(), 2);
        let h3 = graph.create_handle("CA".as_bytes(), 3);

        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h2, h3));

        let p1 = graph.create_path_handle("P1".as_bytes(), false);
        graph.append_step(&p1, h1);
        graph.append_step(&p1, h2);
        graph.append_step(&p1, h3);

        // Generate the forward
        let total_length = find_graph_seq_length(&graph);
        let mut seq_bv: BitVec = BitVec::new_fill(false, total_length + 1);
        let mut node_ref: Vec<NodeRef> = Vec::new();
        let forward = find_sequence_po(&graph, &mut seq_bv, &mut node_ref);

        assert_eq!(total_length, 8);
        assert_eq!("ACGTTTCA", forward);
        assert_eq!(
            *node_ref.get(1).unwrap(),
            NodeRef {
                seq_idx: 3,
                edge_idx: 1,
                edges_to_node: 1
            }
        );
        assert_eq!(
            *node_ref.get(2).unwrap(),
            NodeRef {
                seq_idx: 6,
                edge_idx: 2,
                edges_to_node: 1
            }
        );

        let kmers_on_graph_rust_ver = generate_kmers_linearly(&graph, 3, Some(100), Some(100));
        let kmers_on_graph_dozyg = generate_kmers(&graph, 3, Some(100), Some(100));

        assert_eq!(kmers_on_graph_rust_ver.len(), 10);
        assert_eq!(kmers_on_graph_dozyg.len(), 10);

        for kmer in &kmers_on_graph_rust_ver {
            assert!(kmers_on_graph_dozyg.contains(kmer));
        }
    }

    #[test]
    fn test_self_loop() {
        let mut graph = HashGraph::new();

        //          v--
        // ACG -> TTT | -> CA

        let h1 = graph.create_handle("ACG".as_bytes(), 1);
        let h2 = graph.create_handle("TTT".as_bytes(), 2);
        let h3 = graph.create_handle("CA".as_bytes(), 3);

        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h2, h2));
        graph.create_edge(&Edge(h2, h3));

        let p1 = graph.create_path_handle("P1".as_bytes(), false);
        graph.append_step(&p1, h1);
        graph.append_step(&p1, h2);
        graph.append_step(&p1, h3);

        // Generate the forward
        let total_length = find_graph_seq_length(&graph);
        let mut seq_bv: BitVec = BitVec::new_fill(false, total_length + 1);
        let mut node_ref: Vec<NodeRef> = Vec::new();
        let forward = find_sequence_po(&graph, &mut seq_bv, &mut node_ref);

        assert_eq!(total_length, 8);
        assert_eq!("ACGTTTCA", forward);
        assert_eq!(
            *node_ref.get(1).unwrap(),
            NodeRef {
                seq_idx: 3,
                edge_idx: 1,
                edges_to_node: 2
            }
        );
        // Note that edge_idx is +1 when compared to the previous test, i.e. there is one more edge
        // representing the loop
        assert_eq!(
            *node_ref.get(2).unwrap(),
            NodeRef {
                seq_idx: 6,
                edge_idx: 3,
                edges_to_node: 1
            }
        );

        let kmers_on_graph_rust_ver = generate_kmers_linearly(&graph, 3, Some(100), Some(100));
        let kmers_on_graph_dozyg = generate_kmers(&graph, 3, Some(100), Some(100));

        assert_eq!(kmers_on_graph_rust_ver.len(), 10);
        assert_eq!(kmers_on_graph_dozyg.len(), 10);

        for kmer in &kmers_on_graph_rust_ver {
            assert!(kmers_on_graph_dozyg.contains(kmer));
        }
    }

    #[test]
    fn test_table() {
        let mut graph = create_simple_graph();

        let seq_length = find_graph_seq_length(&graph);

        // Generate the forward
        let total_length = find_graph_seq_length(&graph);
        let mut seq_bv: BitVec = BitVec::new_fill(false, total_length + 1);
        let mut node_ref: Vec<NodeRef> = Vec::new();
        let forward = find_sequence_po(&graph, &mut seq_bv, &mut node_ref);

        let kmers_on_graph = generate_kmers(&graph, 3, Some(100), Some(100));

        let kmers_positions_on_ref: Vec<KmerPos> =
            generate_pos_on_ref(&graph, &kmers_on_graph, &seq_length, &node_ref);

        let hash_build = RandomState::with_seeds(0,0,0,0);
        let hashes = generate_hash(&kmers_on_graph, &hash_build);
        let kmers_mphf = NoKeyBoomHashMap::new_parallel(hashes.clone(), kmers_positions_on_ref.clone());

        for i in 0..kmers_on_graph.len() {
            let kmer = kmers_on_graph.get(i).unwrap();
            let kmer_pos_on_ref = kmers_positions_on_ref.get(i).unwrap();
            let hash = hashes.get(i).unwrap();

            let table_value = kmers_mphf.get(hash).unwrap();
            assert_eq!(kmer_pos_on_ref, table_value);

        }

    }

    #[test]
    fn test_serialization() {
        let mut graph = create_simple_graph();

        let seq_length = find_graph_seq_length(&graph);

        // Generate the forward
        let total_length = find_graph_seq_length(&graph);
        let mut seq_bv: BitVec = BitVec::new_fill(false, total_length + 1);
        let mut node_ref: Vec<NodeRef> = Vec::new();
        let forward = find_sequence_po(&graph, &mut seq_bv, &mut node_ref);

        // Check node_ref serialization
        let encoded_node_ref = bincode::serialize(&node_ref).unwrap();
        let deserialized_node_ref : Vec<NodeRef> = bincode::deserialize(&encoded_node_ref).unwrap();
        assert_eq!(node_ref, deserialized_node_ref);

        // Check bitvec serialization
        let encoded_seq_bv = bincode::serialize(&seq_bv).unwrap();
        let deserialized_seq_bv : BitVec = bincode::deserialize(&encoded_seq_bv).unwrap();
        assert_eq!(seq_bv, deserialized_seq_bv);

        let kmers_on_graph = generate_kmers(&graph, 3, Some(100), Some(100));

        // Check kmer_graph serialization
        let encoded_kmer_graph = bincode::serialize(&kmers_on_graph).unwrap();
        let deserialized_kmer_graph : Vec<Kmer> = bincode::deserialize(&encoded_kmer_graph).unwrap();
        assert_eq!(kmers_on_graph, deserialized_kmer_graph);

        let kmers_positions_on_ref: Vec<KmerPos> =
            generate_pos_on_ref(&graph, &kmers_on_graph, &seq_length, &node_ref);

        let hash_build = RandomState::with_seeds(0,0,0,0);
        let hashes = generate_hash(&kmers_on_graph, &hash_build);
        let kmers_table = NoKeyBoomHashMap::new_parallel(hashes.clone(), kmers_positions_on_ref.clone());

        //Check kmer_table serialization
        let encoded_table = bincode::serialize(&kmers_table).unwrap();
        let deserialized_table : NoKeyBoomHashMap<u64, KmerPos> = bincode::deserialize(&encoded_table).unwrap();
        //assert_eq!(kmers_table, deserialized_table);

        for i in 0..kmers_on_graph.len() {
            let kmer = kmers_on_graph.get(i).unwrap();
            let kmer_pos_on_ref = kmers_positions_on_ref.get(i).unwrap();
            let hash = hashes.get(i).unwrap();

            let table_value = kmers_table.get(hash).unwrap();
            let deserialized_table_value = deserialized_table.get(hash).unwrap();
            assert_eq!(table_value, deserialized_table_value);

        }

    }
}
