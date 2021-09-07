#!/usr/bin/env python
import argparse
import pandas as pd
from glob import glob
from os.path import join, basename, abspath
from re import compile, match

def import_sample_sheet(fp):
    sample_sheet = pd.read_csv(fp,
                               header=None,
                               names=['Sample', 'i7','i5'])

    return(sample_sheet)

def parse_filename(fname, pattern):
    m = match(pattern, fname)

    if m:
        seq_info = {'flowcell': m.groups()[0],
                    'sample': m.groups()[1],
                    'i7': m.groups()[2],
                    'i5': m.groups()[3],
                    'read': m.groups()[4]}
        return(seq_info)
    else:
        return(False)

def parse_seq_files(seq_fps,
                    pattern):
    seq_dict = {}
    unmatched = []
    for file in seq_fps:
        seq_info = parse_filename(basename(file), pattern)

        if seq_info:
            seq_dict[abspath(file)] = seq_info
        else:
            unmatched.append(file)
    return(seq_dict, unmatched)

def match_pe_seqs(sample_df, seq_df, match_col=None):
    seqs = {}
    unmatched = []
    for i, row in sample_df.iterrows():
        if match_col:
            sample = row[match_col]
        else:
            sample = i
        sample_reads = seq_df.loc[seq_df['sample'] == sample,:]
        if len(sample_reads) == 2:
            seqs[sample] = {'r1': sample_reads.loc[sample_reads['read'] == 'R1'].index[0],
                            'r2': sample_reads.loc[sample_reads['read'] == 'R2'].index[0]}
        else:
            unmatched.append(sample)
    return(seqs)

def format_fofn(pe_seqs,
                runtype='paired-end',
                extra=None):
    

    fofn_df = pd.DataFrame.from_dict(pe_seqs,
                                     orient='index')

    fofn_df.index.name = 'sample'

    if extra:
        fofn_df['extra'] = extra
    else:
        fofn_df['extra'] = ''
    fofn_df['runtype'] = runtype
    print(fofn_df)
    fofn_df = fofn_df[['runtype','r1','r2','extra']]
    return(fofn_df)


parser = argparse.ArgumentParser()
parser.add_argument('--sample_sheet', 
                    '-s',
                    required=True,
                    type=str)
parser.add_argument('--match_col', 
                    '-m',
                    required=False,
                    type=str)
parser.add_argument('--output_fp', 
                    '-o',
                    required=True,
                    type=str)
parser.add_argument('--ilm_regex',
                    '-r',
                    type=str,
                    default='^\d+_\d+_\d+_(\w+)_(.+?)_([AGCT]+)_([AGCT]+)_([R12]+)(.+)$')
parser.add_argument('seq_dir',
                    type=str)


def main():
    args = parser.parse_args()

    sample_sheet_fp = args.sample_sheet
    files_dir = args.seq_dir
    ilm_regex = args.ilm_regex
    fofn_fp = args.output_fp
    match_col = args.match_col

    sample_df = import_sample_sheet(sample_sheet_fp)

    files_list = glob(join(files_dir,'*'))

    ilm_re = compile(ilm_regex)
    seq_dict, unmatched = parse_seq_files(files_list, ilm_re)
    seq_df = pd.DataFrame.from_dict(seq_dict, orient='index')
    pe_seqs = match_pe_seqs(sample_df, seq_df, match_col=match_col)
    fofn_df = format_fofn(pe_seqs)

    fofn_df.to_csv(fofn_fp, sep='\t')


if __name__ == "__main__":
    main()

