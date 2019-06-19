#!/usr/bin/env python

import click
import h5py
import numpy as np
import numpy.testing
import pandas as pd
import selene_sdk.sequences

@click.command()
@click.option('-i', '--input_file', 'input_file', type=str, required=True,
                help='the location of the balanced sequences sets for classification')
@click.option('-o', '--output_file', 'output_file', type=str, required=True,
                help='name of the output file for the encoded hdf5 data')
@click.option('-c', '--class_column', 'class_column', type=str, required=True,
                help='the name of the column in the data file with the target class info')
@click.option('-s', '--sequence_column', 'sequence_column', type=str,
                default='sequence', show_default=True,
                help='the name of the column in the data file with sequence data')
def all_seq_save_hdf5(input_file, output_file, class_column, sequence_column):
    """encode DNA sequences and classification column and save to hdf5"""

    print(f'reading cds data from tsv file {input_file}')
    df = pd.read_csv(input_file, sep='\t').filter(items=[sequence_column, class_column])

    print(f'encoding {df.shape[0]} sequences')
    bases_arr = np.array(['A', 'C', 'G', 'T'])
    bases_encoding = { 'A': 0, 'C': 1, 'G': 2, 'T': 3 }
    df[sequence_column] = df[sequence_column].apply(lambda x:
        selene_sdk.sequences.sequence_to_encoding(x, bases_encoding, bases_arr))

    print('creating final data frame with encoded taxonomy flags')
    df['target'] = np.array(df[class_column], dtype=int)
    df.drop(class_column, axis=1, inplace=True)

    print('saving output to {output_file}')
    hdf5 = h5py.File(outfile, 'w')
    sequences = np.array(df.iloc[range(start, end)][sequence_column].values.tolist())
    hdf5.create_dataset('sequence', data=sequences)
    del sequences
    targets = np.array([df.iloc[range(start, end)][class_column].values.tolist()])
    hdf5.create_dataset('target', data=targets)
    hdf5.close()


if __name__ == '__main__':
    all_seq_save_hdf5()
