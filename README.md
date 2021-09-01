# rs-vgaligner
Aligns reads to a Variation Graph by exploiting **Chaining** and **Partial Order Alignment (POA)**.

## How it works
There are three main steps in this program:

1. **Indexing** - We first build an index over the input graph, which contains all of the kmers (of a given size k) encoded by it, and
   their positions w.r.t. the linearized forward/reverse sequence. These positions are stored in a memory-efficient way by 
   using a [Minimal Perfect Hash Function](https://github.com/10XGenomics/rust-boomphf).
2. **Mapping** - We then split the query into kmers, and find perfect maches between these kmers and the index ("**anchors**"). 
   Anchors that are close together are grouped into approximate alignments ("**chains**"), with the constraint that an anchor 
   can appear in at most one chain.
3. **Alignment** - In order to convert chains into **base-level alignments**, we align the subgraph implied by the chain
   to (parts of) the query. We use [abPOA](https://github.com/yangao07/abPOA) (or more specifically [its Rust bindings](https://github.com/HopedWall/rs-abPOA))
   to perform the sequence-to-graph alignment.

## Usage
This crate is mostly inteded to be used as an executable, 
but you can also use its functions if you may want to do so.

### Index

Generate the index with the following command:

```
cargo run --release -- index -k {kmer-size} -i {path-to-gfa-file} -o {output-folder}
```

The index is composed of multiple files, so you might want to add a prefix
after the output folder (i.e. ```-o ./output/my-prefix```). 

There are also additional parameters to limit the number of nodes/edges to be traversed, you can use 
```cargo run -- index --help``` to see them all.

### Mapping and Alignment
TODO
