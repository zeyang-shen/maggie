#!/usr/bin/env python
import sys
sys.path.append('..')

from maggie import score

import numpy as np
import pandas as pd

import os
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate highest motif scores and score differences')
    parser.add_argument("--motifPath",
                        help="path to the motif files",
                        default='/home/zes017/Spacing/Data/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar/', 
                        type=str)
    parser.add_argument("-m", "--motifs", 
                        help="spcify motifs to compute; multiple motifs should be separated by comma like 'SPI1,CEBPB'",
                        type=str)
    parser.add_argument("inputPath", 
                        help="path to a fasta file containing sequences",
                        type = str)
    parser.add_argument("-o", "--output", 
                        help="name of the output file; default has an ending '_motifScores.tsv'",
                        type=str)
    parser.add_argument("-p", "--processor", 
                        help="number of processors to run",
                        default=1,
                        type=int)
    parser.add_argument("--saveDiff", "-D", 
                        help="option to save the differences of highest motif scores", 
                        action="store_true")
    parser.add_argument("--reverse", "-R", 
                        help="the fasta file is organized as negative->positive instead of positive->negaitve as default", 
                        action="store_true")
    args = parser.parse_args()
    
    fasta_file = args.inputPath
    motif_dir = args.motifPath
    score_file = args.output
    proc = args.processor
    saveDiff_bool = args.saveDiff
    rev = args.reverse
    
    if score_file is None:
        score_file = '.'.join(fasta_file.split('.')[:-1])+'_motifScores.tsv'
 
    # Read in motif files
    motif_dict = score.load_motifs(motif_dir)
    if args.motifs:
        motif_list = [i for i in args.motifs.split(',')]
    else:
        motif_list = list((motif_dict.keys()))
    
    # Read in sequences
    id_seq_dict = score.read_fasta(fasta_file)
    sorted_ids = sorted(id_seq_dict.keys())
    
    # Compute motif scores in parallel
    motif_score_df = score.compute_scores_parallel(motif_dict, motif_list, id_seq_dict, proc)
    # Save motif scores
    motif_score_df.to_csv(score_file, sep='\t')
    
    # Compute score difference
    if saveDiff_bool:
        if rev:
            motifScore_diff = score.compute_score_difference(motif_score_df, reverse=True)
        else:
            motifScore_diff = score.compute_score_difference(motif_score_df, reverse=False)
        # Save motif score differences
        diff_file = '.'.join(fasta_file.split('.')[:-1])+'_motifScoresDiff.tsv'
        motifScore_diff.to_csv(diff_file, sep='\t')
    
    