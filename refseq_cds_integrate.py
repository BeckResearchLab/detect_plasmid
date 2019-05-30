#!/usr/bin/env python

import sys

import click
import numpy as np
import pandas as pd
import scipy.io


@click.command()
@click.option('-o', '--output_file', 'output_file', type=str, required=True,
            help='name of the output file with filtered and / or trimmed sequences')
@click.option('-r', '--refseq_input_file', 'refseq_input_file', type=str, required=True,
            help='the location of the refseq taxonomy annotated sequences')
@click.option('-p', '--plasmid_input_file', 'plasmid_input_file', type=str, required=True,
            help='the location of the plasmid taxonomy annotated sequences')
@click.option('-m', '--min_seq_len', 'min_seq_length', type=int, default=0, show_default=True,
            help='the minimum sequence length to be included for integration')
def refseq_cds_integrate(output_file, refseq_input_file, plasmid_input_file, min_seq_length):
    """Integrate refseq and plasmid CDS into a single file"""

    print(f'reading cds tsv data file {plasmid_input_file}')
    df_plasmid = pd.read_csv(plasmid_input_file, sep='\t')

    print(f'reading cds tsv data file {refseq_input_file}')
    df_refseq = pd.read_csv(refseq_input_file, sep='\t')
    df_refseq['is_plasmid'] = 0

    df = df_plasmid.append(df_refseq)

    if min_seq_length > 0:
        print(f'filtering {df.shape[0]} sequences for minimum sequence length of {min_seq_length}')
        df = df.loc[df.sequence.str.len() >= min_seq_length]

    print(f'saving {df.shape[0]} samples to {output_file}')
    df.to_csv(output_file, sep='\t', index=False)


if __name__ == '__main__':
    refseq_cds_integrate()
