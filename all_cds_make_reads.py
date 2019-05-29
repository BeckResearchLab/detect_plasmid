#!/usr/bin/env python

import random
import sys

from Bio import SeqIO
from Bio.Seq import Seq
import click
import numpy as np
import pandas as pd
import scipy.io


@click.command()
@click.option('-r', '--read_length', 'read_length', default=300, type=int,
            show_default=True,
            help='the length of generated reads from the sequences; '+
            'this is accomplished by randomly sampling strings of trim length'+
            'from the full sequence as forward, reverse, complement and '+
            'reverse complement')
@click.option('-m', '--min_seq_length', 'min_seq_length', default=0, type=int,
            show_default=True,
            help='should sequences shortern than this be discarded (0 = disabled)')
@click.option('-o', '--output_file', 'output_file', type=str, required=True,
            help='name of the output file with reads')
@click.option('-i', '--input_file', 'input_file', type=str, required=True,
            help='the location of the taxonomy annotated sequences')
@click.option('-s', '--random_seed', 'random_seed', type=int, default=42,
            show_default=True,
            help='random seed to be used for shuffling and sampling data partitions')
def all_cds_filter(read_length, min_seq_length, output_file, input_file, random_seed):
    """Trim and / or filter sequences by their length"""

    random.seed(random_seed)

    print(f'reading cds tsv data file {input_file}')
    df = pd.read_csv(input_file, sep='\t')

    if min_seq_length == 0:
        min_seq_length = read_length
    print(f'filtering {df.shape[0]} sequences for minimum sequence length of {min_seq_length}')
    df = df.loc[df.sequence.str.len() > min_seq_length]

    df_reads = pd.DataFrame().reindex_like(df)
    print(f'randomly sampling {df.shape[0]} sequences to maximum length of {read_length}')
    df.apply(make_reads, axis=1, max_length=read_length, df_reads=df_reads)

    print(f'saving {df.shape[0]} trimmed and filtered samples to {output_file}')
    df_reads.to_csv(output_file, sep='\t', index=False)


def make_reads(row, max_length, df_reads):
    seq = row['sequence']
    product_id = row['product_id']
    n = int(len(seq) / max_length)
    #print(n, len(seq), max_length, len(seq) - max_length)
    # 4 = forward, reverse, complement, reverse complement
    starts = random.choices(range(len(seq) - max_length), k=n*4)
    for i, start in enumerate(starts):
        read = str(seq[start:start+max_length])
        if i % 4 == 0:
            row['sequence'] = read
        elif i % 4 == 1:
            row['sequence'] = read[::-1]
        elif i % 4 == 2:
            seqr = Seq(read)
            complement = str(seqr.complement())
            row['sequence'] = complement
        else:
            seqr = Seq(read)
            reverse_complement = str(seqr.reverse_complement())
            row['sequence'] = reverse_complement
        row['product_id'] = product_id + '_' + str(i)
        df_reads.append(row)
    return None


if __name__ == '__main__':
    all_cds_filter()
