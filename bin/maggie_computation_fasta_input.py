#!/usr/bin/env python
import sys
sys.path.append('..')

from maggie import score, utils

import numpy as np
import pandas as pd

import os
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Maggie computation pipeline dealing with FASTA files')
    parser.add_argument("originalFile", 
                        help="file1,file2,file3,... fasta file(s) that contain positive sequences; multiple files should be separated by comma without space",
                        type=str)
    parser.add_argument("mutatedFile", 
                        help="file1,file2,file3,... fasta file(s) that contain negative sequences that should be aligned with positive sequences as pairs",
                        type=str)
    parser.add_argument("--motifPath",
                        help="path to the motif files",
                        default='../JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar/', 
                        type=str)
    parser.add_argument("-m", "--motifs", 
                        help="spcify motifs to compute; multiple motifs should be separated by comma like 'SPI1,CEBPB'",
                        type=str)
    parser.add_argument("-o", "--output", 
                        help="name of the output file",
                        default='maggie_output.tsv',
                        type=str)
    parser.add_argument("-p", "--processor", 
                        help="number of processors to run",
                        default=1,
                        type=int)
    args = parser.parse_args()
    
    orig_file = args.originalFile
    mut_file = args.mutatedFile
    motif_dir = args.motifPath
    output_file = args.output
    proc = args.processor
    
    # read in sequences
    if len(orig_file.split(',')) > 1:
        try:
            cat_fasta_files(orig_file.split(','), mut_file.split(','))
            orig_seq_dict = utils.read_fasta('concat_files_ref.fa')
            mut_seq_dict = utils.read_fasta('concat_files_mut.fa')
        except:
            print('Failed to concatenate multiple fasta files.')
    else:
        orig_seq_dict = utils.read_fasta(orig_file)
        mut_seq_dict = utils.read_fasta(mut_file)
        
    if len(orig_seq_dict) != len(mut_seq_dict):
        raise Exception('# positive sequences unequal to # negative sequences. Make sure your inputs come as pairs')
    
    # Read in motif files
    motif_dict = score.load_motifs(motif_dir)
    if args.motifs:
        motif_list = [i for i in args.motifs.split(',')]
    else:
        motif_list = list((motif_dict.keys()))
    
    # Run Maggie pipeline
    results = score.test_all_motifs(motif_dict, orig_seq_dict, mut_seq_dict, top_site=1, p=proc, motif_list=motif_list)
    results.to_csv(output_file, sep='\t')
    
    # Combine similar motifs
    merge_stats = score.combine_similar_motifs(results, 0.6)
    sig_motifs = utils.FDR_cutoff(np.abs(merge_stats['Median p-val']), alpha=0.05)
    merge_stats.loc[sig_motifs].to_csv(output_file.split('.')[0]+'_mergedSignificant.'+output_file.split('.')[1], sep='\t')
    
    
    