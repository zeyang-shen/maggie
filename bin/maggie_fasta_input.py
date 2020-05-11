#!/usr/bin/env python
import sys
sys.path.append('.')

from maggie import score, utils, visual

import numpy as np
import pandas as pd

import os
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='MAGGIE framework working with FASTA files')
    parser.add_argument("posFile", 
                        help="file1,file2,file3,... fasta file(s) that contain positive sequences; multiple files should be separated by comma without space",
                        type=str)
    parser.add_argument("negFile", 
                        help="file1,file2,file3,... fasta file(s) that contain negative sequences that should have the same sequence identifiers as positive sequences to form pairs",
                        type=str)
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
                        type=float)
    parser.add_argument("-sCut", 
                        help="cutoff for calling significance based on FDR values, default = 0.05",
                        default=0.05,
                        type=float)
    parser.add_argument("-T", 
                        help="number of top motif scores to be used to compute for representative motif score, default = 1",
                        default=1,
                        type=int)
    parser.add_argument("--saveDiff", 
                        help="Flag for saving motif score differences. This file can be large, default = False",
                        action='store_true')
    parser.add_argument("--linear", 
                        help="Flag for linear model, default = False",
                        action='store_true')
    parser.add_argument("-p", 
                        help="number of processors to run",
                        default=1,
                        type=int)
    args = parser.parse_args()
    
    orig_file = args.posFile
    mut_file = args.negFile
    motif_dir = args.motifPath
    output_dir = args.output
    mCut = args.mCut
    sCut = args.sCut
    top = args.T
    save = args.saveDiff
    linear = args.linear
    if linear:
        print('Analyzing with linear model')
    proc = args.p
    
    # read in sequences
    if len(orig_file.split(',')) > 1:
        try:
            orig_seq_dict, mut_seq_dict = utils.cat_fasta_files(orig_file.split(','), mut_file.split(','))
        except:
            print('Failed to concatenate multiple fasta files.')
            sys.exit(1)
        print('Successfully concetenated fasta files!')
    else:
        orig_seq_dict = utils.read_fasta(orig_file)
        mut_seq_dict = utils.read_fasta(mut_file)
        
    if len(orig_seq_dict) != len(mut_seq_dict):
        sys.exit('ERROR: unequal # sequences in input files! Make sure your inputs come as pairs')

    
    # Read in motif files
    motif_dict = score.load_motifs(motif_dir)
    if args.motifs:
        motif_list = []
        for i in args.motifs.split(','):
            for k in motif_dict.keys():
                if k.split('$')[0].upper() == i.upper():
                    motif_list.append(k)
        motif_list = np.unique(motif_list)
    else:
        motif_list = list((motif_dict.keys()))
    
    # Create folder to store outputs
    try:
        os.mkdir(output_dir)
    except FileExistsError:
        sys.exit('ERROR: Specified folder exists! Please name a different folder to save results')
    except OSError:
        sys.exit('ERROR: Please check if your specified path exists!')
    
    # Run Maggie pipeline
    print('Running MAGGIE on %d motifs for %d sequences with %d parallel process' % 
          (len(motif_list), len(orig_seq_dict), proc))
    results = score.test_all_motifs(motif_dict, orig_seq_dict, mut_seq_dict, 
                                    top_site=top, p=proc, motif_list=motif_list, linear=linear)
    
    if not linear:
        results.iloc[:,:-1].to_csv(output_dir+'/maggie_output.tsv', sep='\t')
        if save:
            utils.save_raw_data(results, output_dir+'/score_differences.npy')
        # Combine similar motifs
        merge_stats = score.combine_similar_motifs(results, mCut)
        merge_stats.to_csv(output_dir+'/maggie_output_merged.tsv', sep='\t')
        sig_motifs = utils.FDR_cutoff(np.abs(merge_stats['Median p-val']), alpha=sCut)
        merge_stats.loc[sig_motifs].to_csv(output_dir+'/maggie_output_mergedSignificant.tsv', sep='\t')

        # Visualization files
        visual.save_distribution(results, folder=output_dir)
        visual.save_logos(motif_dict, folder=output_dir)
        visual.generate_html(folder=output_dir)
    else:
        results.to_csv(output_dir+'/maggie_output.tsv', sep='\t')
    print('Results are ready in %s' % (output_dir))
    