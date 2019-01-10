import numpy as np
import pandas as pd
import multiprocessing as mp

import os

import Bio
from Bio import motifs
from Bio import SeqIO


def read_fasta(fasta_file):
    '''
    Read in sequences
    '''
    alphabet = Bio.Seq.IUPAC.Alphabet.IUPAC.IUPACUnambiguousDNA() # need to use this alphabet for motif score calculation
    id_seq_dict = {} # {sequenceID: fastq sequence}
    for seq_record in SeqIO.parse(fasta_file, "fasta"):  
        seq_record.seq.alphabet = alphabet
        id_seq_dict[seq_record.id] = seq_record.seq
        
    return id_seq_dict


def load_motifs(motif_dir="/home/zes017/Spacing/Data/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar/"):
    '''read in motifs: motifs have to be in jaspar format'''
    motif_dict = {}
    nuc = ['A', 'C', 'G', 'T']
    for mf in os.listdir(motif_dir):
        with open(motif_dir + mf) as f:
            m = motifs.read(f, 'jaspar')
            counts = np.array([m.counts[n] for n in nuc])
            m.pseudocounts = np.mean(counts.sum(axis=0))//10
            m.background = None
            motif_dict[m.name] = m
            
    return motif_dict


def compute_scores(motif_dict, motif, seq_dict):
    fwd_pssm = motif_dict[motif].pssm
    rev_pssm = fwd_pssm.reverse_complement()
    scores = []
    sorted_ids = sorted(seq_dict.keys())
    
    for sid in sorted_ids:
        seq = seq_dict[sid]
        fwd_scores = fwd_pssm.calculate(seq) # scores for forward orientation
        rev_scores = rev_pssm.calculate(seq) # scores for reverse orientation
            
        # get the highest score in forward direction
        try:
            max_fwd_score = np.max(fwd_scores)
        except ValueError:
            max_fwd_score = fwd_pssm.min
        # get the highest score in reverse direction
        try:
            max_rev_score = np.max(rev_scores)
        except ValueError:
            max_rev_score = rev_pssm.min

        # determine which orientation is better
        if max_fwd_score > max_rev_score:
            scores.append(max_fwd_score)
        else:
            scores.append(max_rev_score)
            
    return motif, scores


def compute_scores_parallel(motif_dict, motif_list, seq_dict, p=1):
    '''
    Compute motif scores in parallel
    '''
    pool = mp.Pool(processes=p)
    results = [pool.apply_async(compute_scores, args=(motif_dict, motif, seq_dict)) for motif in motif_list]
    results = [p.get() for p in results]
    sorted_ids = sorted(seq_dict.keys())
    
    # Combine results from multiple processors
    motif_score_dict = {}
    for r in results:
        motif, score = r
        motif_score_dict[motif] = score
    motif_score_df = pd.DataFrame(motif_score_dict, index=sorted_ids)
    
    return motif_score_df


def compute_score_difference(motifScore_df, reverse=False):
    group_region = motifScore_df.T.groupby(np.arange(len(motifScore_df))//2, axis=1)
    if reverse:
        motifScore_diff = group_region.diff().T.iloc[1::2]
    else:
        motifScore_diff = -group_region.diff().T.iloc[1::2]
    drop_bools = np.any(motifScore_diff.isnull(), axis=1)
    print('# drop:', sum(drop_bools))
    motifScore_diff = motifScore_diff.loc[~drop_bools]
    
    return motifScore_diff


def merge_score_diff(score_diff_files, cut_percent=0):
    print('Merging files:')
    all_files_df = pd.DataFrame()
    for file in score_diff_files:
        print('--', file)
        diff_df = pd.read_csv(file, sep='\t', index_col=0)
        all_files_df = pd.concat([all_files_df, diff_df.T], axis=1)
    motifs = all_files_df.index.values
    
    # filter out distrurbing outliers of each motif
    print('Filtering out', cut_percent*100, '% outliers')
    X = []
    for i in range(len(all_files_df)):
        df = all_files_df.iloc[i]
        cut_df = pd.qcut(df, [cut_percent, 1-cut_percent])
        X.append(np.array(df.loc[cut_df.notnull()]))
    X = np.array(X)
    
    return X, motifs
    
    