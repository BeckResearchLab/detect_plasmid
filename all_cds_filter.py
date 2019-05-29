#!/usr/bin/env python

import random
import sys

import click
import numpy as np
import pandas as pd
import scipy.io


@click.command()
@click.option('-t', '--trim_seq_length', 'trim_length', default=0, type=int,
            show_default=True,
            help='should sequences be trimmed down to a length (0 = disabled)'+
            'this is accomplished by randomly sampling strings of trim length'+
            'from the full sequence')
@click.option('-m', '--min_seq_length', 'min_seq_length', default=0, type=int,
            show_default=True,
            help='should sequences shortern than this be discarded (0 = disabled)')
@click.option('-o', '--output_file', 'output_file', type=str, required=True,
            help='name of the output file with filtered and / or trimmed sequences')
@click.option('-i', '--input_file', 'input_file', type=str, required=True,
            help='the location of the taxonomy annotated sequences')
@click.option('-r', '--random_seed', 'random_seed', type=int, default=42,
            show_default=True,
            help='random seed to be used for shuffling and sampling data partitions')
def all_cds_filter(trim_length, min_seq_length, output_file, input_file, random_seed):
    """Trim and / or filter sequences by their length"""

    random.set(random_seed)

    if trim_length <= 0 and min_seq_length <= 0:
        print('no trim length or minimum sequence filter length specified')
        print('exiting because there is nothing to do!')
        sys.exit(0)

    print(f'reading cds tsv data file {input_file}')
    df = pd.read_csv(input_file, sep='\t')
    if min_seq_length > 0: 
        print(f'filtering {df.shape[0]} sequences for minimum sequence length of {min_seq_length}')
        df = df.loc[df.sequence.str.len() >= min_seq_length]
    if trim_length > 0:
        print(f'randomly sampling {df.shape[0]} sequences to maximum length of {trim_length}')
        df["sequence"] = df["sequence"].apply(random_substring_sample, trim_length)

    print(f'saving {df.shape[0]} trimmed and filtered samples to {output_file}')
    df.to_csv(output_file, sep='\t', index=False)


def random_substring_sample(seq, max_length):
    i = random.sample(range(len(seq) - max_length), 1)
    return seq[i:i+max_length]


if __name__ == '__main__':
    all_cds_filter()
