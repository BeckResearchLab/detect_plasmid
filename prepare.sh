#!/bin/bash

PLASMIDS_PATH=/work/data/NCBI_plasmids/plasmid
REFSEQ_PATH=/work/data/refseq
THREADS=24
MIN_SEQ_LEN=1000
MAX_PLASMID_LEN=500000
MIN_REFSEQ_LEN=1000000
RANDOM_SEED=42
# no limit on samples, use all positives
CLASS_SAMPLES=0

if [ ! -e plasmid_seq.tsv ]; then
	echo "extracting plasmid sequences"
	./plasmid_seq_extractor.py --plasmids_fna $PLASMIDS_PATH/ncbi_plasmid.fna \
			--plasmids_gbff $PLASMIDS_PATH/ncbi_plasmid.gbff \
			--output_file plasmid_seq.tsv --max_length $MAX_PLASMID_LEN
fi

if [ ! -e refseq_seq.tsv ]; then
	echo "extracting refseq sequences"
	./refseq_seq_extractor.py --refseq_path $REFSEQ_PATH \
			--output_file refseq_seq.tsv --threads $THREADS \
			--min_length $MIN_REFSEQ_LEN
fi

if [ ! -e all_seq.tsv ]; then
	echo "integrating refseq and plasmid seqs while length filtering"
	./all_seq_integrate.py --refseq_input_file refseq_seq.tsv \
			--plasmid_input_file plasmid_seq.tsv --min_seq_len $MIN_SEQ_LEN \
			--output_file all_seq.tsv
fi

if [ ! -e all_reads.tsv ]; then
	# this samples sequences like reads, but not necessarily as reads
	./all_seq_make_reads.py --read_length $MIN_SEQ_LEN --input_file all_seq.tsv \
			--output_file all_reads.tsv --random_seed $RANDOM_SEED
fi

if [ ! -e balanced_reads.tsv ]; then
	echo "balancing representation of classes"
	./all_seq_balance.py --input_file all_reads.tsv --output_file balanced_reads.tsv \
			--positive_samples $CLASS_SAMPLES --random_seed $RANDOM_SEED
fi

if [ ! -e all_seq_train.h5 -o ! -e all_seq_valid.h5 -o ! -e all_seq_test.h5 ]; then
	echo "encoding and saving data to hdf5 files"
	./all_seq_save_hdf5.py --input_file balanced_reads.tsv \
			--train_frac 0.91 --valid_frac 0.01 --test_frac 0.08 \
			--train_file all_seq_train.h5 --valid_file all_seq_valid.h5 \
			--test_file all_seq_test.h5
fi
