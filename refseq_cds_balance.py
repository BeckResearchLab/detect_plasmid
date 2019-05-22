#!/usr/bin/env python

import click
import numpy as np
import pandas as pd


@click.command()
@click.option('-l', '--tax_level', 'tax_level', required=True,
        type=click.Choice(['kingdom', 'phylum', 'class', 'order',
                            'family', 'genus']),
        help='what taxonomy level should the data be partitioned on')
@click.option('-t', '--taxa', 'taxa', type=str, required=True,
        help='the taxonomy class that will be balanced for sampling')
@click.option('-o', '--output_file', 'output_file', type=str, required=True,
        help='name of the output file containing balanced samples')
@click.option('-i', '--input_file', 'input_file', type=str, required=True,
        help='the location of the taxonomy annotated sequences')
@click.option('-r', '--random_seed', 'random_seed', type=int, default=42,
        show_default=True,
        help='random seed to be used for shuffling and sampling data partitions')
@click.option('-s', '--positive_samples', 'positive_samples', type=int,
        default=0, show_default=True,
        help='limit the number of positive samples to this number (disabled = 0)')
def refseq_cds_balance(tax_level, taxa, input_file, output_file, random_seed, positive_samples):
    """Balance samples of sequences for a binary classification at a given taxonomy level and class"""

    print(f'reading cds tsv data file {input_file}')
    df = pd.read_csv(input_file, sep='\t').filter(items=['sequence', tax_level])
    print(f'there were {df.shape[0]} sequences in {input_file}')

    print(f'finding sequences at the {tax_level} matching {taxa}')
    positives = df.loc[df[tax_level] == taxa]
    print(f'found {positives.shape[0]} samples')
    if positive_samples > 0:
        print(f'down sampling positive samples to {positive_samples} while shuffling')
        positives = positives.sample(n=positive_samples, random_state=random_seed)
    print(f'finding sequences not matching {tax_level} if {taxa}')
    negatives = df.loc[df[tax_level] != taxa]
    print(f'found {negatives.shape[0]} samples to be randomly sampled down to {positives.shape[0]} samples')
    negatives = negatives.sample(n=positives.shape[0], random_state=random_seed)
    print(f'concatenating positive and negative samples')
    df = positives.append(negatives)
    # attempt to allow garbage collection before shuffling
    positives = None
    negatives = None
    print(f'shuffling the order of samples')
    df = df.sample(frac=1, random_state=random_seed).reset_index(drop=True)
    print(f'saving {df.shape[0]} balanced samples to {output_file}')
    df.to_csv(output_file, sep='\t', index=False)


if __name__ == '__main__':
    refseq_cds_balance()
