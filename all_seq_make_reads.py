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
            help='should sequences shorter than this be discarded (0 = disabled)')
@click.option('-o', '--output_file', 'output_file', type=str, required=True,
            help='name of the output file with reads')
@click.option('-i', '--input_file', 'input_file', type=str, required=True,
            help='the location of the taxonomy annotated sequences')
@click.option('-s', '--random_seed', 'random_seed', type=int, default=42,
            show_default=True,
            help='random seed to be used for shuffling and sampling data partitions')
@click.option('-c', '--coverage', 'coverage', type=int, default=1,
            show_default=True,
            help='coverage like parameter that is used to boost sampling of sequences')
@click.option('-n', '--neg_coverage', 'neg_coverage', type=int, default=0,
            show_default=True,
            help='coverage like parameter that is used to boost sampling of negative class sequences')
@click.option('-p', '--pos_coverage', 'pos_coverage', type=int, default=0,
            show_default=True,
            help='coverage like parameter that is used to boost sampling of positive class sequences')
@click.option('-a', '--class_column', 'class_column', type=str, required=True,
            show_default=True,
            help='column name to use for differential coverage / sampling of class sequences')
def all_seq_make_reads(read_length, min_seq_length, output_file, input_file, random_seed, coverage, pos_coverage, neg_coverage, class_column):
    """Trim and / or filter sequences by their length"""

    random.seed(random_seed)

    print(f'reading seq tsv data file {input_file}')
    df = pd.read_csv(input_file, sep='\t')

    if min_seq_length == 0:
        min_seq_length = read_length
    print(f'filtering {df.shape[0]} sequences for minimum sequence length of {min_seq_length}')
    df = df.loc[df.sequence.str.len() > min_seq_length]
    update_product_id = False
    if 'product_id' in df.columns:
        update_product_id = True

    print(f'generating reads from {df.shape[0]} sequences with length of {read_length}')
    if coverage != 1:
        print(f'default coverage is {coverage}x')
    if pos_coverage != 0:
        print(f'positive sample coverage is {pos_coverage}x')
    if neg_coverage != 0:
        print(f'negative sample coverage is {neg_coverage}x')
    if pos_coverage !=0 or neg_coverage != 0:
        print(f'class column has the name {class_column}')
    df_reads = make_reads(df, read_length, update_product_id, coverage, pos_coverage, neg_coverage, class_column)

    print(f'saving {df_reads.shape[0]} reads to {output_file}')
    df_reads.drop(columns=['index'])
    df_reads.to_csv(output_file, sep='\t', index=False)


def make_reads(df_source, max_length, update_product_id, coverage, pos_coverage, neg_coverage, class_column):
    df_reads = pd.DataFrame().reindex_like(df_source)
    df_reads.drop(df_reads.index, inplace=True)
    reads = []
    # not using apply here because in pandas .24 this wouldn't properly reduce
    count = 0
    for index, row in df_source.iterrows():
        seq = row['sequence']
        if update_product_id:
            product_id = str(row['product_id'])
        n = int(len(seq) / max_length)
        # 4 = forward, reverse, complement, reverse complement
        cov = 1
        if coverage == 0:
            cov = coverage
            if pos_coverage != 0 and row[class_column] == 1:
                cov = pos_coverage
            elif neg_coverage != 0 and row[class_column] == 0:
                cov = neg_coverage

        starts = random.choices(range(len(seq) - max_length), k=n*4*cov)
        for i, start in enumerate(starts):
            read_row = row
            if update_product_id:
                read_row['product_id'] = product_id + '_' + str(i)
            read = str(seq[start:start+max_length])
            if i % 4 == 0:
                # forward
                read_row['sequence'] = read
            elif i % 4 == 1:
                # reverse
                read_row['sequence'] = read[::-1]
            elif i % 4 == 2:
                # complement
                seqr = Seq(read)
                complement = str(seqr.complement())
                read_row['sequence'] = complement
            else:
                # reverse complement
                seqr = Seq(read)
                reverse_complement = str(seqr.reverse_complement())
                read_row['sequence'] = reverse_complement
            reads.append(read_row)

            count = count + 1
            if count > 0 and count % 1000000 == 0:
                print(f"{count} reads generated")
    print(f"assembling data frame from {count} reads")
    return df_reads.append(pd.DataFrame(reads, columns=df_source.columns)).reset_index()


if __name__ == '__main__':
    all_seq_make_reads()
