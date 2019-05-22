#!/usr/bin/env python

import sys

import click
import numpy as np
import pandas as pd
import scipy.io


@click.command()
@click.option('-t', '--trim_seq_length', 'trim_length', default=0, type=int,
            show_default=True,
            help='should sequences be trimmed down to a length (0 = disabled)')
@click.option('-m', '--min_seq_length', 'min_seq_length', default=0, type=int,
            show_default=True,
            help='should sequences shortern than this be discarded (0 = disabled)')
@click.option('-o', '--output_file', 'output_file', type=str, required=True,
            help='name of the output file with filtered and / or trimmed sequences')
@click.option('-i', '--input_file', 'input_file', type=str, required=True,
            help='the location of the taxonomy annotated sequences')
def refseq_cds_filter(trim_length, min_seq_length, output_file, input_file):
    """Trim and / or filter sequences by their length"""

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
        print(f'trimming {df.shape[0]} sequences to maximum length of {trim_length}')
        df["sequence"] = df.sequence.str.slice(stop=trim_length)

    print(f'saving {df.shape[0]} trimmed and filtered samples to {output_file}')
    df.to_csv(output_file, sep='\t', index=False)


if __name__ == '__main__':
    refseq_cds_filter()
