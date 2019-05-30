#!/bin/bash

REFSEQ_PATH=/work/data/refseq
THREADS=24
MIN_SEQ_LEN=300
#CLASS_SAMPLES=100000

./plasmid_cds_extractor.py -f plasmids/ncbi_plasmid.fna -g plasmids/ncbi_plasmid.gbff -o plasmid_cds.tsv
./refseq_cds_integrate.py -r /work/dacb/detect_hgt/300bp/refseq_cds_filtered.tsv -p plasmid_cds.tsv -o all_cds.tsv
./cds_balance
./all_cds_filter

exit



if [ ! -e refseq_cds.tsv ]; then
	./refseq_cds_extractor.py --refseq_path $REFSEQ_PATH \
			--output_file refseq_cds.tsv --threads $THREADS --include_all True
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
