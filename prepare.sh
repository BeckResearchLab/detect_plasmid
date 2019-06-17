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

if [ ! -e balanced_seq.tsv ]; then
	echo "balancing representation of classes"
	./all_seq_balance.py --input_file all_seq.tsv --output_file balanced_seq.tsv \
			--positive_samples $CLASS_SAMPLES --random_seed $RANDOM_SEED
fi

#### SPLIT OUT THE SETS
echo fix me

# training
if [ ! -e training_reads.tsv ]; then
	echo "generating read like sequences for training data"
	# this samples sequences like reads, but not necessarily as reads
	./all_seq_make_reads.py --read_length $MIN_SEQ_LEN 
			--input_file training_seq.tsv \
			--output_file training_reads.tsv --random_seed $RANDOM_SEED \
			--pos_coverage 10 --neg_coverage 1 --class_column 'is_plasmid'
fi

if [ ! -e training.h5 ]; then
	echo "encoding and saving training data to hdf5 files"
	./all_seq_save_hdf5.py --input_file training_reads.tsv \
			--train_frac 1 --valid_frac 0 --test_frac 0 \
			--train_file training.h5 --valid_file /tmp/a.h5 \
			--test_file /tmp/b.h5
fi

# validation
if [ ! -e validation_reads.tsv ]; then
	echo "generating read like sequences for validation data"
	# this samples sequences like reads, but not necessarily as reads
	./all_seq_make_reads.py --read_length $MIN_SEQ_LEN 
			--input_file validation_seq.tsv \
			--output_file validation_reads.tsv --random_seed $RANDOM_SEED \
			--pos_coverage 10 --neg_coverage 1 --class_column 'is_plasmid'
fi

if [ ! -e validation.h5 ]; then
	echo "encoding and saving validation data to hdf5 files"
	./all_seq_save_hdf5.py --input_file validation_reads.tsv \
			--train_frac 0 --valid_frac 1 --test_frac 0 \
			--train_file /tmp/a.h5 --valid_file /tmp/validation.h5 \
			--test_file /tmp/b.h5
fi

# testing
if [ ! -e testing_reads.tsv ]; then
	echo "generating read like sequences for testing data"
	# this samples sequences like reads, but not necessarily as reads
	./all_seq_make_reads.py --read_length $MIN_SEQ_LEN 
			--input_file testing_seq.tsv \
			--output_file testing_reads.tsv --random_seed $RANDOM_SEED \
			--pos_coverage 10 --neg_coverage 1 --class_column 'is_plasmid'
fi

if [ ! -e testing.h5 ]; then
	echo "encoding and saving testing data to hdf5 files"
	./all_seq_save_hdf5.py --input_file testing_reads.tsv \
			--train_frac 0 --valid_frac 0 --test_frac 0 \
			--train_file /tmp/a.h5 --valid_file /tmp/b.h5 \
			--test_file testing.h5
fi
##### 
if [ ! -e balanced_reads.tsv ]; then
	echo "generating read like sequences"
	# this samples sequences like reads, but not necessarily as reads
	./all_seq_make_reads.py --read_length $MIN_SEQ_LEN --input_file all_seq.tsv \
			--output_file balanced_reads.tsv --random_seed $RANDOM_SEED \
			--pos_coverage 10 --neg_coverage 1 --class_column 'is_plasmid'
fi

if [ ! -e all_seq_train.h5 -o ! -e all_seq_valid.h5 -o ! -e all_seq_test.h5 ]; then
	echo "encoding and saving data to hdf5 files"
	./all_seq_save_hdf5.py --input_file balanced_reads.tsv \
			--train_frac 0.91 --valid_frac 0.01 --test_frac 0.08 \
			--train_file all_seq_train.h5 --valid_file all_seq_valid.h5 \
			--test_file all_seq_test.h5
fi
