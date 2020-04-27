#!/usr/bin/env python
import sys
sys.path.append('.')

from maggie import score, utils, visual

import numpy as np
import pandas as pd

from Bio import Seq

import os
import argparse
import requests

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Maggie framework working with VCF file; Installation of HOMER is required!')
    parser.add_argument("vcfFile", 
                        help="VCF file that contains testing variants",
                        type=str)
    parser.add_argument("genome", 
                        help="reference genome",
                        type=str)
    parser.add_argument("-e", "--effect", 
                        help="name of the column in vcfFile for effect sizes that compare alternative vs. reference alleles. If not specified, assume a2 always associated with a higher signal than a1",
                        type=str)
    parser.add_argument("-a1",
                        help="name of the column in vcfFile for testing reference alleles, default = 'REF'",
                        default="REF",
                        type=str)
    parser.add_argument("-a2",
                        help="name of the column in vcfFile for testing alternative alleles, default = 'ALT'",
                        default="ALT",
                        type=str)
    parser.add_argument("-S", "--size",
                        help="size of sequences to test around variants, default = 100",
                        default=100,
                        type=int)
    parser.add_argument("--motifPath",
                        help="path to the motif files",
                        default="./data/JASPAR2020_CORE_vertebrates_motifs/",
                        type=str)
    parser.add_argument("-m", "--motifs", 
                        help="spcify motifs to compute; multiple motifs should be separated by comma without space like 'SPI1,CEBPB'",
                        type=str)
    parser.add_argument("-o", "--output", 
                        help="output directory; by default, will create a new folder under current path",
                        default='./maggie_output/',
                        type=str)
    parser.add_argument("-mCut", 
                        help="cutoff for merging similar motifs; a float value ranging from 0 (merge everything) to 1 (no merging at all), default = 0.6",
                        default=0.6,
                        type=int)
    parser.add_argument("-sCut", 
                        help="cutoff for calling significance based on FDR values, default = 0.05",
                        default=0.05,
                        type=int)
    parser.add_argument("-T", 
                        help="number of top motif scores to be used to compute for representative motif score, default = 1",
                        default=1,
                        type=int)
    parser.add_argument("--save", 
                        help="Flag indicating whether to save the motif score differences for individual motifs. This file can be large, default = False",
                        action='store_true')
    parser.add_argument("-p", 
                        help="number of processors to run",
                        default=1,
                        type=int)
    args = parser.parse_args()
    
    vcf_file = args.vcfFile
    genome = args.genome
    ref_col = args.a1.upper()
    alt_col = args.a2.upper()
    effect_col = args.effect.upper()
    size = args.size
    motif_dir = args.motifPath
    output_dir = args.output
    mCut = args.mCut
    sCut = args.sCut
    top = args.T
    save = args.save
    if not save:
        print('"--save" flag off: will not save score differences for individual motifs')
    proc = args.p
    
    # read in VCF file
    with open(vcf_file, 'r') as rf:
        skips = 0
        for line in rf:
            if line[0] == '#' and len(line.split('\t')) == 1:
                skips += 1
            else:
                break
    vcf_df = pd.read_csv(vcf_file, sep='\t', skiprows=skips)
    vcf_df.columns = [str(col).upper() for col in vcf_df.columns.values]
    if ref_col not in vcf_df.columns.values or alt_col not in vcf_df.columns.values:
        sys.exit('ERROR: reference and/or alternative alleles are not defined!')
    if effect_col not in vcf_df.columns.values:
        sys.exit('ERROR: cannot find effect sizes! Please make sure you have the column name correct')

    # step 1: get surrounding sequences
    if len(genome.split('.')) == 1: # download reference genome
        genome_path = './data/genomes/'+genome+'.fa'
        if not os.path.exists(genome_path):
            print('Downloading reference genome', genome)
            try:
                url = 'http://homer.ucsd.edu/zeyang/maggie/genomes/'+genome+'.fa'
                r = requests.get(url, stream=True)
                r.raise_for_status()
                total_size = int(r.headers.get('content-length', 0))
                block_size = 1024
                try:
                    from tqdm import tqdm # load package to show the progress bar
                    t=tqdm(total=total_size, unit='B', unit_scale=True)
                    with open(genome_path, 'wb') as f:
                        for data in r.iter_content(block_size):
                                t.update(len(data))
                                f.write(data)
                        t.close()
                except:
                    with open(genome_path, 'wb') as f:
                        f.write(r.content)
            except:
                print('ERROR: failed to download the specified genome')
                raise
    else:
        genome_path = genome
    
    # read reference genome
    try:
        print('Reading reference genome:', os.path.abspath(genome_path))
        genomes = utils.load_genome(genome_path)
    except FileExistsError:
        print('ERROR: please check if your specified genome path exists!')
        raise
    
    pos_data = utils.data_prep(vcf_file, genomes, size=size, skiprows=skips+1, file_format='vcf')
    chr_list = ['chr'+str(ch) if 'chr' not in str(ch) else str(ch) for ch in vcf_df.iloc[:,0]]
    mids = vcf_df.iloc[:,1]
    seq_keys = [c+'.'+str(mids[k]) for k,c in enumerate(chr_list)]
    seq_dict = dict()
    for i in range(len(pos_data)):
        seq_dict[seq_keys[i]] = pos_data[i][0]
    
    # step 2: generate positive and negative sequences
    vcf_df.index = seq_keys
    alphabet = Seq.IUPAC.Alphabet.IUPAC.IUPACUnambiguousDNA()
    high_seqs = []
    low_seqs = []
    indices = []
    for snp in seq_dict.keys():
        ref_allele = vcf_df.loc[snp][ref_col].upper()
        alt_allele = vcf_df.loc[snp][alt_col].upper()
        if str(ref_allele) == 'nan':
            print('Skipping:', snp)
            continue

        seq = str(seq_dict[snp]).upper()
        if ref_allele == seq[len(seq)//2-1:len(seq)//2-1+len(ref_allele)]:
            ref = seq
            alt = seq[:len(seq)//2-1]+alt_allele+seq[len(seq)//2-1+len(ref_allele):]
        elif alt_allele == seq[len(seq)//2-1:len(seq)//2-1+len(alt_allele)]:
            ref = seq[:len(seq)//2-1]+ref_allele+seq[len(seq)//2-1+len(alt_allele):]
            alt = seq
        else:
            print('Skipping:', snp)
            continue

        indices.append(snp)
        ref = Seq.Seq(ref, alphabet=alphabet)
        alt = Seq.Seq(alt, alphabet=alphabet)
        t_val = vcf_df.loc[snp][effect_col]
        # save high vs. low allele
        if t_val > 0:
            high_seqs.append(alt)
            low_seqs.append(ref)
        else:
            high_seqs.append(ref)
            low_seqs.append(alt)
    seq_df = pd.DataFrame(index=indices, columns=['high_seq', 'low_seq'])
    seq_df['high_seq'] = high_seqs
    seq_df['low_seq'] = low_seqs
    orig_seq_dict = dict(seq_df['high_seq'])
    mut_seq_dict = dict(seq_df['low_seq'])
    
    # step 3: run MAGGIE
    # Read in motif files
    motif_dict = score.load_motifs(motif_dir)
    if args.motifs:
        motif_list = [i for i in args.motifs.split(',')]
    else:
        motif_list = list((motif_dict.keys()))
    
    # Create folder to store outputs
    try:
        os.mkdir(output_dir)
    except FileExistsError:
        pass
    except OSError:
        print('Please check if your specified output path exists!')
        raise
    
    # run pipeline
    print('Running MAGGIE on %d motifs for %d sequences with %d parallel process' % 
          (len(motif_list), len(orig_seq_dict), proc))
    results = score.test_all_motifs(motif_dict, orig_seq_dict, mut_seq_dict, top_site=top, p=proc, motif_list=motif_list)
    if save:
        results.to_csv(output_dir+'/maggie_output.tsv', sep='\t')
    else:
        results.iloc[:,:-1].to_csv(output_dir+'/maggie_output.tsv', sep='\t')
    
    # Combine similar motifs
    merge_stats = score.combine_similar_motifs(results, mCut)
    merge_stats.to_csv(output_dir+'/maggie_output_merged.tsv', sep='\t')
    sig_motifs = utils.FDR_cutoff(np.abs(merge_stats['Median p-val']), alpha=sCut)
    merge_stats.loc[sig_motifs].to_csv(output_dir+'/maggie_output_mergedSignificant.tsv', sep='\t')
    
    # Visualization files
    visual.save_distribution(results, folder=output_dir)
    visual.save_logos(motif_dict, folder=output_dir)
    visual.generate_html(folder=output_dir)
    print('Results are ready in %s' % (output_dir))
    