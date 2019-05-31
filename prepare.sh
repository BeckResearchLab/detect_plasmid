#!/bin/bash

REFSEQ_PATH=/work/data/refseq
THREADS=24
MIN_SEQ_LEN=300
RANDOM_SEED=42
# no limit on samples, use all positives
CLASS_SAMPLES=0


if [ ! -e plasmid_cds.tsv ]; then
	./plasmid_cds_extractor.py --plasmids_fna plasmids/ncbi_plasmid.fna --plasmids_gbff plasmids/ncbi_plasmid.gbff \
			--output_file plasmid_cds.tsv
fi

if [ ! -e all_cds.tsv ]; then
	./refseq_cds_integrate.py --refseq_input_file /work/dacb/detect_hgt/refseq_cds.tsv \
			--plasmid_input_file plasmid_cds.tsv --min_seq_len $MIN_SEQ_LEN \
			--output_file all_cds.tsv
fi

if [ ! -e balanced_cds.tsv ]; then
	./all_cds_balance.py --input_file all_cds.tsv --output_file balanced_cds.tsv \
			--positive_samples $CLASS_SAMPLES --random_seed $RANDOM_SEED
fi

if [ ! -e balanced_reads.tsv ]; then
	./all_cds_make_reads.py --read_length $MIN_SEQ_LEN --input_file balanced_cds.tsv \
			--output_file balanced_reads.tsv --random_seed $RANDOM_SEED
fi

if [ ! -e all_cds_train.h5 -o ! -e all_cds_valid.h5 -o ! -e all_cds_test.h5 ]; then
	./all_cds_save_hdf5.py --input_file balanced_reads.tsv \
			--train_frac 0.91 --valid_frac 0.01 --test_frac 0.08 \
			--train_file all_cds_train.h5 --valid_file all_cds_valid.h5 \
			--test_file all_cds_test.h5
fi
