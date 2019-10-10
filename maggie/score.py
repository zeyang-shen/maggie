import os
import multiprocessing as mp

import numpy as np
import pandas as pd

import Bio
from Bio import motifs
from Bio import SeqIO

from scipy.stats import ttest_1samp, wilcoxon


def load_motifs(motif_dir, pseudocounts=0.01):
    '''
    read in motifs; motifs have to be in jaspar format as below:
    
        >MA0002.2       RUNX1
        A  [   287    234    123     57      0     87      0     17     10    131    500 ]
        C  [   496    485   1072      0     75    127      0     42    400    463    158 ]
        G  [   696    467    149      7   1872     70   1987   1848    251     81    289 ]
        T  [   521    814    656   1936     53   1716     13     93   1339   1325   1053 ]
    
    Parameters:
        motif_dir: folder that contains motif files; one file for individual motif
        pseudocounts: w.r.t. position weight matrix, the probability adding to every nucleotide
    '''
    motif_dict = {}
    nuc = ['A', 'C', 'G', 'T']
    for mf in os.listdir(motif_dir):
        with open(motif_dir + mf) as f:
            m = motifs.read(f, 'jaspar')
            counts = np.array([m.counts[n] for n in nuc])
            avg_counts = counts.sum(axis=0).mean()
            m.pseudocounts = avg_counts*pseudocounts
            m.background = None
            motif_dict[m.name] = m
            
    return motif_dict


def compute_scores(bio_motif, seq_dict, top_site):
    '''
    compute motif scores across sequences and 
    output top scores to represent log-likelihood of being bound by transcription factor
    
    Parameters:
        bio_motif: Bio motif object used to compute motif scores
        seq_dict: python dictionary containing sequences in IUPAC alphabet
        top_site: the number of top motif scores output for downstream analysis
    
    Outputs:
        (motif name, motif scores, positions of scores on sequence)
    '''
    fwd_pssm = bio_motif.pssm
    rev_pssm = fwd_pssm.reverse_complement()
    scores = []
    pos = []
    sorted_ids = sorted(seq_dict.keys())
    
    for sid in sorted_ids:
        seq = seq_dict[sid]
        fwd_scores = fwd_pssm.calculate(seq) # scores for forward orientation
        rev_scores = rev_pssm.calculate(seq) # scores for reverse orientation
        
        concat_scores = list(fwd_scores) + list(rev_scores)
        if len(concat_scores) < 0:
            max_scores = [0]
            max_pos = [-1]
        else:
            max_pos = np.argsort(concat_scores)[::-1][:top_site]
            max_scores = np.array(concat_scores)[max_pos]
        max_scores = np.nan_to_num(max_scores, 0)
        max_scores = max_scores*(max_scores > 0)
        scores.append(np.log2(np.sum(2**max_scores))) #additive log-likelihood
#         scores.append(np.sum(max_scores))
        pos.append(max_pos%len(fwd_scores))
            
    return bio_motif.name, scores, pos


def test_one_motif(bio_motif, orig_seq_dict, mut_seq_dict, top_site):
    '''
    test for score differences of one motif
    positive and negative sequences are required to have aligned IDs
    
    Parameters:
        bio_motif: Bio motif object used to compute motif scores
        orig_seq_dict: python dictionary containing positive sequences
        mut_seq_dict: python dictionary containing negative sequences
        top_site: the number of top motif scores output for downstream analysis
    
    Outputs:
        (motif name, statistic, p value, score differences)
    '''
    orig_score = compute_scores(bio_motif, orig_seq_dict, top_site)
    mut_score = compute_scores(bio_motif, mut_seq_dict, top_site)
    score_diff = list(np.array(orig_score[1]) - np.array(mut_score[1]))
    pv_ = []
    stat_ = []
    for k in range(1000):
        score_diff_sampled = np.random.choice(score_diff, replace=True, size=len(score_diff))
        tmp_stat, tmp_pv = wilcoxon(score_diff_sampled)
        pv_.append(-np.log10(tmp_pv)*np.sign(np.median(score_diff_sampled[score_diff_sampled!=0])))
        stat_.append(tmp_stat)
    return (bio_motif.name, stat_, pv_, score_diff)


def test_all_motifs(motif_dict, orig_seq_dict, mut_seq_dict, top_site=1, p=1, motif_list=None):
    '''
    test for all motifs in the given dictionary
    positive and negative sequences are required to have aligned IDs
    
    Parameters:
        motif_dict: Bio motif object used to compute motif scores
        orig_seq_dict: python dictionary containing positive sequences
        mut_seq_dict: python dictionary containing negative sequences
        top_site: the number of top motif scores output for downstream analysis
        p: specify th number of cores to do parrelel processes
    
    Outputs:
        A pandas dataframe including all the results
    '''
    if not motif_list:
        motif_list = np.sort(list(motif_dict.keys()))
    pool = mp.Pool(processes=p)
    results = [pool.apply_async(test_one_motif, args=(motif_dict[m], orig_seq_dict, mut_seq_dict, top_site)) 
               for m in motif_list]
    results = [r.get() for r in results]
    # format results to panda dataframe
    results_df = pd.DataFrame(results)
    results_df.columns = ['motif', 'stats list', 'p-val list', 'score difference']
    results_df = results_df.set_index('motif')
    results_df['5% stats'] = [np.percentile(r, 5) for r in results_df['stats list']]
    results_df['Median stats'] = [np.percentile(r, 50) for r in results_df['stats list']]
    results_df['95% stats'] = [np.percentile(r, 95) for r in results_df['stats list']]
    results_df['5% p-val'] = [np.percentile(r, 5) for r in results_df['p-val list']]
    results_df['Median p-val'] = [np.percentile(r, 50) for r in results_df['p-val list']]
    results_df['95% p-val'] = [np.percentile(r, 95) for r in results_df['p-val list']]
    results_df = results_df[['5% stats', 'Median stats', '95% stats', 
                             '5% p-val', 'Median p-val', '95% p-val', 
                             'score difference', 'stats list', 'p-val list']]
    results_df = results_df.fillna(0)
    
    return results_df


def merge_tuple(lsts):
    '''
    Merge tuples
    '''
    sets = [set(lst) for lst in lsts if lst]
    merged = True
    while merged:
        merged = False
        results = []
        while sets:
            common, rest = sets[0], sets[1:]
            sets = []
            for x in rest:
                if x.isdisjoint(common):
                    sets.append(x)
                else:
                    merged = True
                    common |= x
            results.append(common)
        sets = results
    return sets


def combine_similar_motifs(input_df, similarity_cutoff=0.6):
    '''
    merge the results of similar motifs to reduce the long list
    
    Input:
        input_df: the results dataframe output from function "test_all_motifs"
    
    Parameters:
        similarity_cutoff: cutoff correlation coefficient for calling similar motifs (default: 0.6)
    
    Outputs:
        Average statistics and p-values for merged sets of similar motifs
    '''
    motif_names = input_df.index.values
    score_diffs = np.array([r for r in input_df['score difference']])
    corr_df = pd.DataFrame(np.corrcoef(score_diffs), index=motif_names, columns=motif_names)
    
    tuple_list = np.array(corr_df[corr_df < 1][corr_df > similarity_cutoff].stack().index)
    merge_sets = merge_tuple(tuple_list)
    rest_list = [{s} for s in set(motif_names)-set([element for tup in tuple_list for element in tup])]
    merge_sets += rest_list
    
    merge_stats = pd.DataFrame()
    for s in merge_sets:
        mean_set = input_df.loc[list(s)].iloc[:,:6].mean(axis=0)
        mean_set.name = '|'.join(list(s))
        merge_stats = pd.concat([merge_stats, mean_set], axis=1, sort=True)
    merge_stats = merge_stats.T

    return merge_stats    
    