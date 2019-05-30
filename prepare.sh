#!/bin/bash

REFSEQ_PATH=/work/data/refseq
THREADS=24
MIN_SEQ_LEN=300
CLASS_SAMPLES=100000



if [ ! -e plasmid_cds.tsv ]; then
	./plasmid_cds_extractor.py --plasmids_fna plasmids/ncbi_plasmid.fna --plasmids_gbff plasmids/ncbi_plasmid.gbff \
			--output_file plasmid_cds.tsv
fi

if [ ! -e all_cds.tsv ]; then
	./refseq_cds_integrate.py --refseq_input_file /work/dacb/detect_hgt/${MIN_SEQ_LEN}bp/refseq_cds_filtered.tsv \
			--plasmid_input_file plasmids_cds.tsv --min_seq_len 300 \
			--output_file all_cds.tsv
fi

if [ ! -e balanced_cds.tsv ]; then
	./all_cds_balance.py --input_file all_cds.tsv --output_file balanced_cds.tsv \
			--positive_samples $CLASS_SAMPLES --random_seed 42
fi

exit

if [ ! -e refseq_cds_train.mat -o ! -e refseq_cds_valid.mat -o ! -e refseq_cds_test.mat ]; then
	./refseq_cds_savemat.py --input_file refseq_cds_filtered_balanced.tsv \
			--tax_level family --taxa Enterobacteriaceae \
			--train_frac 0.7 --valid_frac 0.2 --test_frac 0.1 \
			--train_file refseq_cds_train.mat --valid_file refseq_cds_valid.mat \
			--test_file refseq_cds_test.mat
fi
