# rs-vgaligner
Aligns reads to a Variation Graph by combining **Seed-and-Extend approaches** with **Partial Order Alignment (POA)**.

## How it works
There are three main steps in this program:

1. **Indexing** - We first build an index over the input graph, which contains all of the kmers (of a given size k) encoded by it, and
   their positions w.r.t. the linearized forward/reverse sequence. These positions are stored in a memory-efficient [Minimal Perfect Hash Function](https://github.com/10XGenomics/rust-boomphf).
2. **Mapping** - We then split the query into kmers, and find perfect maches between these kmers and the linearization ("**anchors**"). 
   Anchors that are close together are grouped into approximate alignments ("**chains**"). This step is inspired by [Minimap2](https://github.com/lh3/minimap2).
3. **Alignment** - In order to convert chains into **base-level alignments**, we align the subgraph implied by the chain
   to the query. We use [abPOA](https://github.com/yangao07/abPOA) (or more specifically [its Rust bindings](https://github.com/HopedWall/rs-abPOA))
   to perform the sequence-to-graph alignment.

## Usage
This crate is mostly inteded to be used as an executable, but you can also use its functions if you may want to do so.

In order to install vgaligner as an executable, clone this repository, move to its main folder and run:

```
cargo install --path .
```

Graphs must be sorted in order to ensure that vgaligner works correctly. You can do so with [odgi](https://github.com/pangenome/odgi):

```
odgi sort -i unsorted_graph.gfa -o - -p Ygs -P | odgi view -i - -g > sorted_graph.gfa
```

### Index

Generate the index with the following command:

```
vgaligner index -k {kmer-size} -i {input.gfa} -o {index.idx}
```

There are also additional parameters to limit the number of nodes/edges to be traversed, you can use 
```vgaligner index --help``` to see them all.

### Mapping and Alignment

Map reads to the graph with the following command:

```
vgaligner map -i {index.idx} -f {reads.fa/fq} -o {output.gaf}
```

This will output the chains in GAF format. 

If you also want to perform 
the alignment with abPOA (this is highly recommended!), you will need to pass the ```--also-align``` parameter. 
This will generate an additional GAF file which contains the alignment.

You can see all of the available parameters with ```vgaligner map --help```.

## Validation
This tool has been tested on graphs obtained from [HLA-zoo](https://github.com/ekg/HLA-zoo). 

The validation pipeline uses **snakemake**, and can be executed by using the following commands:

```
cd experiments-snakemake && snakemake -k --cores {cores_to_be_used}
```

We recommend having at least 16 GB RAM and a CPU with 4+ cores.
