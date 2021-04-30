# rs-vgaligner
Aligns reads to a Variation Graph by exploiting **Partial Order Alignment (POA)**.

## How it works
There are three main steps in this program:

1. **Indexing**:
   we first build an index over the input graph, which contains all of the kmers (of a given size k) encoded by it, and
   their positions w.r.t. the linearized forward/reverse sequence. These positions are stored in a memory-efficient way by 
   using a Minimal Perfect Hash Function (MPHF).
2. **Mapping** - TODO
3. **Querying** - TODO