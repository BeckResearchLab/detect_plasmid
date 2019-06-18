#!/usr/bin/env python

import random
import sys

import click
import numpy as np
import pandas as pd
import scipy.io


@click.command()
@click.option('-f', '--min_fragment_length', 'min_fragment_length', type=int,
            required=True,
            help='when breaking up sequences larger than this length, this is the minimum size for a fragment')
@click.option('-o', '--output_file', 'output_file', type=str, required=True,
            help='name of the output file with filtered and / or trimmed sequences')
@click.option('-i', '--input_file', 'input_file', type=str, required=True,
            help='the location of the taxonomy annotated sequences')
@click.option('-r', '--random_seed', 'random_seed', type=int, default=42,
            show_default=True,
            help='random seed to be used for shuffling and sampling data partitions')
@click.option('-s', '--sequence_column', 'sequence_column', type=str, default='sequence',
            show_default=True,
            help='name of the column containing the sequence data')
@click.option('-c', '--class_column', 'class_column', type=str, default='is_plasmid',
            show_default=True,
            help='name of the column containing the class information')
def all_seq_fragment(output_file, input_file, random_seed, min_fragment_length, sequence_column, class_column):
    """Fragment sequences into chunks with a minimum length"""

    random.seed(random_seed)

    print(f'reading cds tsv data file {input_file}')
    df = pd.read_csv(input_file, sep='\t')
    print(f'chopping sequences longer than {min_fragment_length}bp into')
    print(f'    sequencings of at least {min_fragment_length} bp long')
    print(f'    saving the output fo {output_file}')
    f = open(output_file, 'w')
    f.write(f'{sequence_column}\t{class_column}\n')

    fragments = 0
    for index, row in df.iterrows():
        is_class = row[class_column]
        sequence = row[sequence_column]
        # fragment
        if len(sequence) < min_fragment_length:
            f.write(f'{sequence}\t{is_class}\n')
            fragments = fragments + 1
            continue
        for fragment in fragment_sequence(sequence, min_fragment_length):
            f.write(f'{fragment}\t{is_class}\n')
            fragments = fragments + 1

    f.close()
    print(f'saved {fragments} fragments')


def fragment_sequence(sequence, min_fragment_length):
    lseq = len(sequence)
    num_chunks = lseq // min_fragment_length
    chunk_size = lseq // num_chunks

    if lseq % num_chunks:
        chunk_size += 1

    iterator = iter(sequence)
    for _ in range(num_chunks):
        accumulator = list()
        for _ in range(chunk_size):
            try: accumulator.append(next(iterator))
            except StopIteration: break
        yield ''.join(accumulator)


if __name__ == '__main__':
    all_seq_fragment()
