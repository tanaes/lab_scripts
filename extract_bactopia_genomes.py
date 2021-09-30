#!/usr/bin/env python

import click
import pandas as pd
from os import makedirs
from os.path import join
from shutil import copyfile


@click.command()
@click.option('-s', '--summary', 'summary_fp',
              required=True,
              help='path Bactopia summary file')
@click.option('-b', '--bactopia', 'bactopia_dir',
              required=True,
              help='path to Bactopia base dir')
@click.option('-o', '--output', 'out_dir',
              required=True,
              help='path to output directory')
@click.option('-e', '--copy_exclude', 'copy_exclude',
              required=False,
              default=False,
              help='also copy "exclude" genomes')
def extract_genomes(summary_fp, bactopia_dir, out_dir, copy_exclude):
    summary_df = pd.read_csv(summary_fp, sep='\t', header=0, index_col=0)

    assem_fp_pt = join(bactopia_dir, '%s/assembly/%s.fna')
    makedirs(out_dir, exists_ok=True)

    for i, row in summary_df.iter_rows():
        rank = row['rank']

        if rank != 'exclude' or copy_exclude:
            copyfile(assem_fp_pt % i, join(out_dir, '%s.fna' % i))


if __name__ == '__main__':
    extract_genomes()
