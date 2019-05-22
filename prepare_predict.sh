#!/bin/bash

mkdir -f predict_sequences
head -1000 refseq_cds_filtered.tsv | awk -F'\t' '{ printf(">%s\t%s\n%s\n", $9, $7, $10); }' > predict_sequences/seqs.fa
