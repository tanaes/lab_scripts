#!/usr/bin/env python

import click
import gzip
import pandas as pd
import numpy as np
from copy import deepcopy
from os.path import join, splitext
from os import makedirs
import pysam
from glob import glob
from os.path import basename


def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

def load_genome(genome_fp):
    genome = {}
    with open(genome_fp, 'r') as f:
        for name, seq, _ in readfq(f):
            genome[name] = list(seq)
    return(genome)

def parse_genome_variants(genome_snp_df, genome, min_freq=0.8, min_cov=5):

    genome_cons = deepcopy(genome)
    genome_mvar = deepcopy(genome)
    genome_nvar = deepcopy(genome)

    cons_changes = 0
    var_changes = 0
    n_changes = 0
    
    mut_types = {}

    for i, row in genome_snp_df.iterrows():
        scaffold = row['scaffold']
        pos = row['position']
        
        mut = row['mutation_type']
        if mut not in mut_types:
            mut_types[mut] = 1
        else:
            mut_types[mut] += 1
        
        if row['con_freq'] == 1:
            genome_cons[scaffold][pos] = row['con_base']
            genome_mvar[scaffold][pos] = row['con_base']
            genome_nvar[scaffold][pos] = row['con_base']
            cons_changes += 1
        elif row['con_freq'] >= min_freq:
            genome_cons[scaffold][pos] = 'N'
            genome_mvar[scaffold][pos] = row['con_base']
            genome_nvar[scaffold][pos] = row['var_base']
            var_changes += 1
        else:
            genome_cons[scaffold][pos] = 'N'
            genome_mvar[scaffold][pos] = 'N'
            genome_nvar[scaffold][pos] = 'N'
            n_changes += 1

    
    return(genome_cons, genome_mvar, genome_nvar, mut_types)

def write_genome_fasta(fp,
                       genome):
    with open(fp, 'w') as f:
        for contig in genome:
            seq = ''.join(genome[contig])
            f.write('>{0}\n{1}\n'.format(contig, seq))

def split_genome(host,
                 genome_name,
                 genome_snp_df,
                 genome, 
                 outdir,
                 min_freq = 0.8,
                 min_cov = 5,):
    """
    returns:
        mut_types
    """
    genome_cons, genome_mvar, genome_nvar, mut_types =  parse_genome_variants(genome_snp_df, genome, min_freq=0.8)
    
    return(genome_cons, genome_mvar, genome_nvar, mut_types)
    
def genome_cov(genome, bam_fp, min_cov=5):
    bam = pysam.AlignmentFile(bam_fp, "rb")
    for contig in genome:
        cov = [x.n for x in bam.pileup(contig)]
        for i, n in enumerate(cov):
            if n < min_cov:
                genome[contig][i] = 'N'

def split_sample_genomes(genome_df,
                         snp_df,
                         genome_fps,
                         bam_fp,
                         outdir,
                         mg_ind,
                         min_cov=5,
                         min_freq=0.8):

    mut_info = {}

    for genome_name in genome_fps.index:
        if genome_name not in genome_df.index:
            print('Genome {0} not in outuput\n'.format(genome_name))
            continue
        if genome_df.loc[genome_name, 'coverage'] < 5:
            print('Genome {0} coverage insufficient\n'.format(genome_name))
            continue
        genome_snp_df = snp_df.loc[snp_df['genome'] == genome_name]
        genome = load_genome(genome_fps[genome_name])
        ref_sp = genome_df.loc[genome_name, 'host_species']
        ref_ind = genome_df.loc[genome_name, 'host_individual']
        
        print('splitting genome: %s' % genome_name)
        genome_cons, genome_mvar, genome_nvar, mut_types = split_genome(mg_ind,
                                 genome_name,
                                 genome_snp_df,
                                 genome, 
                                 outdir,
                                 min_freq=min_freq,
                                 min_cov=min_cov)
        
        parsed = {'cons': genome_cons,
                  'maj': genome_mvar,
                  'min': genome_nvar}
        
        print('applying coverage threshold')
        for g in parsed:
            genome_cov(parsed[g],
                       bam_fp,
                       min_cov=min_cov)
            write_genome_fasta(join(outdir, '{0}.{2}.{1}{3}'.format(mg_ind,
                                                                    g,
                                                                    *splitext(genome_name))),
                               parsed[g])
    
        mut_info[genome_name] = mut_types

    genome_data_df = pd.concat([genome_df, pd.DataFrame.from_dict(mut_info,orient='index')],
                               axis=1)
    
    return(genome_data_df)


@click.command()
@click.option('-a',
              '--align_dir',
              'align_dir',
              required=True,
              help='path to inStrain alignment directory')
@click.option('-g',
              '--genome_list_fp',
              'genome_list_fp',
              required=True,
              help='path to file of genome fps')
@click.option('-s',
              '--stb_fp',
              'stb_fp',
              required=True,
              help='path to stb file')
@click.option('-p',
              '--profile_dir',
              'profile_dir',
              required=True,
              help='path to inStrain profiles directory')
@click.option('-o',
              '--output_dir',
              'output_dir',
              required=True,
              help='path to output directory')
def extract_genomes(align_dir,
                    profile_dir,
                    genome_list_fp,
                    stb_fp,
                    output_dir):

    makedirs(output_dir, exist_ok=True)

    mgs = glob(join(profile_dir,'*'))
    genome_fps = pd.read_csv(genome_list_fp,
                             sep='\t',
                             header=None,
                             index_col=0,
                             squeeze=True)


    stb_df = pd.read_csv(stb_fp, 
                         sep='\t',
                         header=None,
                         names=['scaffold','genome'],
                         index_col=0)

    for mg in mgs:
        print('looking at %s' % mg)
        
        mg_name = basename(mg)
        mg_sp = mg_name.split('_')[0]
        mg_ind = mg_name.split('.')[0]

        genome_df_fp = join(mg, 'output', mg_name + '_genome_info.tsv')
        genome_df = pd.read_csv(genome_df_fp, header=0, index_col=0, sep='\t')
        
        genome_df['host_species'] = [x.split('-')[1].split('_')[0] for x in genome_df.index]
        genome_df['host_individual'] = ['_'.join(x.split('-')[1].split('_')[:2]) for x in genome_df.index]

        bam_fp = join(align_dir, '{0}.sorted.bam'.format(mg_ind))
        
        snp_fp = join(mg, 'output', '{0}_SNVs.tsv.gz'.format(mg_name))
        snp_df = pd.read_csv(snp_fp, compression='gzip', sep='\t')
        snp_df = pd.merge(left=snp_df, right=stb_df, left_on='scaffold', right_index=True)

        genome_data_df = split_sample_genomes(genome_df,
                                              snp_df,
                                              genome_fps,
                                              bam_fp,
                                              output_dir,
                                              mg_ind)

        print(mg_ind)
        genome_data_df['mg_host_individual'] = mg_ind
        genome_data_df['mg_host_species'] = mg_sp

        genome_data_dfs[mg_ind] = genome_data_df
        
        
    df_list = []
    for mg in genome_data_dfs:
        df = genome_data_dfs[mg]
        df['file'] = df.index
        df_list.append(df)
    total_df = pd.concat(df_list)

    total_df['within_individual'] = False
    total_df.loc[(total_df['mg_host_individual'] == total_df['host_individual']), 
                 'within_individual'] = True
    total_df['within_species'] = False
    total_df.loc[(total_df['mg_host_species'] == total_df['host_species']), 
                 'within_species'] = True

    total_df['comparison'] = 'between_species'
    total_df.loc[total_df['within_species'], 'comparison'] = 'within_species'
    total_df.loc[total_df['within_individual'], 'comparison'] = 'within_individual'

    total_df['dNdS'] = total_df['N']/total_df['S']

    total_df.set_index(['mg_host_individual','host_individual'],
                       inplace=True)
    total_df.to_csv(join(output_dir, 'snp_info.csv'))


if __name__ == '__main__':
    extract_genomes()

