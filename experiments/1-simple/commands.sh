#!/bin/bash 

# Install vgaligner (just for safety)
cargo install --path ../..

# Generate the reads (no headers though)
vg sim -x graph.gfa -n 1 -s 77 > reads.fa

# Add header to reads.fa
awk -v n=1 '{print ">read" n "\n" $0}{read+=1}' reads.fa > reads_header.fa

# Generate the alignments (binary -- gam)
vg sim -x graph.gfa -n 1 -s 77 -a > aln.gam

# Run vgaligner (how do I decide k?)
vgaligner index  -i graph.gfa -k 11
vgaligner map -i graph.idx -f reads_header.fa --also-align

# Convert my result to gam
vg convert graph.gfa -F reads_header-alignments.gaf > reads_header-alignments.gam

# Compare alignment
vg gamcompare reads_header-alignments.gam aln_header.gam

# (optional) Convert the alignment from vg sim to gaf
#vg convert graph.gfa -G aln.gam > aln.gaf

# Add the correct read header to the gaf
#awk -v n=1 '$1 = "read" n' aln.gaf > aln_header.gaf

# Conver to gam
#vg convert graph.gfa -F aln_header.gaf > aln_header.gam

# Compare alignment
#vg gamcompare reads_header-alignments.gam aln_header.gam
