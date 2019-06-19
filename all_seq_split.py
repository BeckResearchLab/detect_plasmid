#!/usr/bin/env python

import click
import numpy as np
import numpy.testing
import pandas as pd
import selene_sdk.sequences

def df_named_subset_save(title, outfile, frac, df, start, end):
    print(f"saving {title} data to {outfile} ({frac * 100.}% = {len(range(start, end))} samples)")
    sdf = df.iloc[range(start, end)]
    sdf.to_csv(outfile, sep='\t', index=False)


@click.command()
@click.option('-i', '--input_file', 'input_file', type=str, required=True,
                help='the location of the balanced sequences sets for classification')
@click.option('-r', '--train_file', 'train_file', type=str, required=True,
                help='name of the output file for the training data set')
@click.option('-v', '--valid_file', 'valid_file', type=str, required=True,
                help='name of the output file for the validation data set')
@click.option('-s', '--test_file', 'test_file', type=str, required=True,
                help='name of the output file for the test data set')
@click.option('-n', '--train_frac', 'train_frac', type=float, default=0.7,
                show_default=True,
                help='the fraction of samples to use for the training data set')
@click.option('-a', '--valid_frac', 'valid_frac', type=float, default=0.2,
                show_default=True,
                help='the fraction of samples to use for the validation data set')
@click.option('-e', '--test_frac', 'test_frac', type=float, default=0.1,
                show_default=True,
                help='the fraction of samples to use for the test data set')
@click.option('-c', '--class_column', 'class_column', type=str, default='is_plasmid',
                show_default=True,
                help='the name of the column in the data file with the target class info')
@click.option('-q', '--sequence_column', 'sequence_column', type=str,
                default='sequence', show_default=True,
                help='the name of the column in the data file with sequence data')
def all_seq_split(input_file, train_file, valid_file, test_file, train_frac, valid_frac, test_frac, class_column, sequence_column):
    """split a set of annoted sequences into training, validation and test sets"""

    numpy.testing.assert_almost_equal(train_frac + valid_frac + test_frac, 1.,
        err_msg='the fractions of training, validation and test data do not equal 1')

    print(f'reading cds data from tsv file {input_file}')
    df = pd.read_csv(input_file, sep='\t').filter(items=['sequence', 'is_plasmid'])
    print(f'randomizing the order of the rows in the dataframe before split')
    df = df.sample(frac=1).reset_index(drop=True)

    print('creating final data frame with just sequence and target columns')
    df[class_column] = np.array(df[class_column], dtype=int)

    print('splitting in training, validation, and test sets')
    max_train = int(train_frac * df.shape[0])
    df_named_subset_save('training', train_file, train_frac, df, 0, max_train)
    valid_i = int(valid_frac * df.shape[0])
    df_named_subset_save('validation', valid_file, valid_frac, df, max_train, max_train+valid_i)
    df_named_subset_save('test', test_file, test_frac, df, max_train+valid_i, df.shape[0])


if __name__ == '__main__':
    all_seq_split()
