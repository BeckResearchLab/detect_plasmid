#!/bin/bash

REFSEQ_PATH=/work/data/refseq
THREADS=24
MIN_SEQ_LEN=750
CLASS_SAMPLES=100000

if [ ! -e refseq_cds.tsv ]; then
	./refseq_cds_extractor.py --refseq_path $REFSEQ_PATH \
			--output_file refseq_cds.tsv --threads $THREADS
fi

if [ ! -e refseq_cds_filtered.tsv ]; then
	./refseq_cds_filter.py --input_file refseq_cds.tsv --output_file refseq_cds_filtered.tsv \
			--min_seq_length $MIN_SEQ_LEN --trim_seq_length $MIN_SEQ_LEN
fi

if [ ! -e refseq_cds_filtered_balanced.tsv ]; then
	./refseq_cds_balance.py --tax_level family --taxa Enterobacteriaceae \
			--output_file refseq_cds_filtered_balanced.tsv \
			--input_file refseq_cds_filtered.tsv \
			--positive_samples $CLASS_SAMPLES --random_seed 42
fi

if [ ! -e refseq_cds_train.mat -o ! -e refseq_cds_valid.mat -o ! -e refseq_cds_test.mat ]; then
	./refseq_cds_savemat.py --input_file refseq_cds_filtered_balanced.tsv \
			--tax_level family --taxa Enterobacteriaceae \
			--train_frac 0.7 --valid_frac 0.2 --test_frac 0.1 \
			--train_file refseq_cds_train.mat --valid_file refseq_cds_valid.mat \
			--test_file refseq_cds_test.mat
fi
