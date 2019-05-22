#!/usr/bin/env python

from datetime import datetime
import io
import os
import shutil

from Bio import SeqIO
import click
import multiprocessing
import numpy as np
import pandas as pd


def process_pool_init(lock, outfile, include_all):
    global output_file_lock
    output_file_lock = lock
    global output_f
    output_f = outfile
    global include_all_loci
    include_all_loci = include_all

def gff_cds_extract(filepath):
    output = io.StringIO()
    try:
        for seq_record in SeqIO.parse(filepath, 'genbank'):
            taxonomy = seq_record.annotations['taxonomy']
            if not include_all_loci and \
                    ('plasmid' in seq_record.description or
                        'extrachromosomal' in seq_record.description):
                continue
            for feature in seq_record.features:
                if feature.type == 'CDS':
                    output.write(f"{filepath}\t{seq_record.id}\t{taxonomy[0] if len(taxonomy) > 0 else np.nan}\t{taxonomy[1] if len(taxonomy) > 1 else np.nan}\t{taxonomy[2] if len(taxonomy) > 2 else np.nan}\t{taxonomy[3] if len(taxonomy) > 3 else np.nan}\t{taxonomy[4] if len(taxonomy) > 4 else np.nan}\t{taxonomy[5] if len(taxonomy) > 5 else np.nan}\t{feature.qualifiers['protein_id'][0] if 'protein_id' in feature.qualifiers else np.nan}\t{feature.location.extract(seq_record).seq}\n")
    except AttributeError:
        print(f'parsing of file {filepath} failed')

    output.seek(0)
    output_file_lock.acquire()
    shutil.copyfileobj(output, output_f)
    output_f.flush()
    output_file_lock.release()


@click.command()
@click.option('-t', '--threads', 'threads', default=16, type=int,
        help='number of parallel Genbank parser threads')
@click.option('-r', '--refseq_path', 'refseq_path', type=str, required=True,
        help='path to the root of the refseq download')
@click.option('-o', '--output_file', 'output_file', type=str, required=True,
        help='name of the output file w/ taxonomy annotations and sequences')
@click.option('-g', '--genbank_postfix', 'genbank_postfix', type=str,
        default='.gbff', show_default=True,
        help='postfix for Genbank files')
@click.option('-i', '--include_all', 'include_all', type=bool,
        default=False, show_default=True,
        help='should plasmids and extrachromosomal elements be included')
def refseq_cds_extractor(threads, refseq_path, output_file, genbank_postfix, include_all):
    """Extract the taxonomy and CDS sequences from a collection of GBFF files"""

    start_time = datetime.now()

    print(f'finding Genbank files in {refseq_path}')
    gbff_files = []
    for root, dirs, files in os.walk(refseq_path):
        for file_ in files:
            filepath = os.path.join(root, file_)
            if filepath.endswith(genbank_postfix):
                gbff_files.append(filepath)

    print(f'scanning {len(gbff_files)} files to extract taxonomy and CDS sequences')
    print(f'using {threads} parallel parsers')

    f = open(output_file, 'w')
    f.write('gff_file\tid\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tproduct_id\tsequence\n')
    f.flush()

    lock = multiprocessing.Lock()
    process_pool = multiprocessing.Pool(threads,
            initializer=process_pool_init, initargs=(lock, f, include_all, ))
    process_pool.map(gff_cds_extract, gbff_files)
    process_pool.close()
    process_pool.join()

    f.close()

    stop_time = datetime.now()
    total_time = stop_time - start_time
    print(f'run time was: {total_time}')


if __name__ == '__main__':
    refseq_cds_extractor()
