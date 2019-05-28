#!/usr/bin/env python

from datetime import datetime
import io
import os
import shutil

from Bio import SeqIO
import click
import numpy as np
import pandas as pd


@click.command()
@click.option('-f', '--plasmids_fna', 'plasmids_fna_path', type=str, required=True,
        help='path to the FASTA nucleotide file of the plasmids')
@click.option('-g', '--plasmids_gbff', 'plasmids_gbff_path', type=str,
        help='path to the Genbank file of the plasmids')
@click.option('-o', '--output_file', 'output_file', type=str, required=True,
        help='name of the output file w/ taxonomy annotations and sequences')
def refseq_cds_extractor(plasmids_fna_path, plasmids_gbff_path, output_file):
    """Extract the taxonomy and CDS sequences from a GBFF and FNA file for plasmids"""

    start_time = datetime.now()

    print(f'reading FASTA file {plasmids_fna_path}')
    print(f'reading GBFF file {plasmids_gbff_path}')

    f = open(output_file, 'w')
    f.write('gff_file\tid\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tproduct_id\tsequence\n')
    f.flush()

    fna_dict = SeqIO.to_dict(SeqIO.parse(plasmids_fna_path, 'fasta'))

    try:
        for seq_record in SeqIO.parse(plasmids_gbff_path, 'genbank'):
            taxonomy = seq_record.annotations['taxonomy']
            for feature in seq_record.features:
                if feature.type == 'CDS':
                    if 'contig' in seq_record.annotations:
                        seq = 'NN'
                        output.write(f"{filepath}\t{seq_record.id}\t{taxonomy[0] if len(taxonomy) > 0 else np.nan}\t{taxonomy[1] if len(taxonomy) > 1 else np.nan}\t{taxonomy[2] if len(taxonomy) > 2 else np.nan}\t{taxonomy[3] if len(taxonomy) > 3 else np.nan}\t{taxonomy[4] if len(taxonomy) > 4 else np.nan}\t{taxonomy[5] if len(taxonomy) > 5 else np.nan}\t{feature.qualifiers['protein_id'][0] if 'protein_id' in feature.qualifiers else np.nan}\t{seq}\n")
                    else:
                        output.write(f"{filepath}\t{seq_record.id}\t{taxonomy[0] if len(taxonomy) > 0 else np.nan}\t{taxonomy[1] if len(taxonomy) > 1 else np.nan}\t{taxonomy[2] if len(taxonomy) > 2 else np.nan}\t{taxonomy[3] if len(taxonomy) > 3 else np.nan}\t{taxonomy[4] if len(taxonomy) > 4 else np.nan}\t{taxonomy[5] if len(taxonomy) > 5 else np.nan}\t{feature.qualifiers['protein_id'][0] if 'protein_id' in feature.qualifiers else np.nan}\t{feature.location.extract(seq_record).seq}\n")
    except AttributeError:
        print(f'parsing of file {filepath} failed')

    f.close()

    stop_time = datetime.now()
    total_time = stop_time - start_time
    print(f'run time was: {total_time}')


if __name__ == '__main__':
    refseq_cds_extractor()
