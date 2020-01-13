import os
import multiprocessing as mp

import numpy as np
import pandas as pd

from Bio import motifs, SeqIO

from scipy.stats import ttest_1samp, wilcoxon


def load_motifs(motif_dir, pseudocounts=0.05):
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
            motif_dict[m.name+'$'+m.matrix_id] = m
            
    return motif_dict


def compute_scores(bio_motif, seq_dict, top_site=1):
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
#         pos.append(max_pos%len(fwd_scores))
            
    return scores


def test_one_motif(bio_motif, orig_seq_dict, mut_seq_dict, top_site=1):
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
    score_diff = np.array(orig_score) - np.array(mut_score)
    # Statistical testing on the original score differences
    nonzero_diff = score_diff[score_diff!=0]
    pv_ = []
    stat_ = []
    if len(nonzero_diff) < 10:
        tmp_stat = 0
        tmp_pv = 1
    else:
        tmp_stat, tmp_pv = wilcoxon(score_diff)
    pv_.append(-np.log10(tmp_pv)*np.sign(np.median(nonzero_diff)))
    stat_.append(tmp_stat)
    # Bootstrapping
    for k in range(999):
        score_diff_sampled = np.random.choice(score_diff, replace=True, size=len(score_diff))
        nonzero_diff = score_diff_sampled[score_diff_sampled!=0]
        if len(nonzero_diff) < 10:
            tmp_stat = 0
            tmp_pv = 1
        else:
            tmp_stat, tmp_pv = wilcoxon(score_diff_sampled)
        pv_.append(-np.log10(tmp_pv)*np.sign(np.median(nonzero_diff)))
        stat_.append(tmp_stat)
    return (bio_motif.name, bio_motif.matrix_id, stat_, pv_, list(score_diff))


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
    # parallel processing to compute motif score differences
    pool = mp.Pool(processes=p)
    try:
        from tqdm import tqdm # load package to show the progress bar
        pbar = tqdm(total=len(motif_list))
        results = [pool.apply_async(test_one_motif, args=(motif_dict[m], orig_seq_dict, mut_seq_dict, top_site), 
                                    callback=lambda _: pbar.update(1)) 
                   for m in motif_list]
    except:
        print('Missing package to show the progress')
        results = [pool.apply_async(test_one_motif, args=(motif_dict[m], orig_seq_dict, mut_seq_dict, top_site)) 
                   for m in motif_list]
    results = [r.get() for r in results]
    pool.close()
    pool.join()
    # format results to panda dataframe
    results_df = pd.DataFrame(results)
    results_df.columns = ['motif', 'id', 'stats list', 'p-val list', 'score difference']
    results_df = results_df.set_index('id')
    results_df['5% stats'] = [np.around(np.percentile(r, 5), decimals=1) for r in results_df['stats list']]
    results_df['Median stats'] = [np.around(r[0], decimals=1) for r in results_df['stats list']]
    results_df['95% stats'] = [np.around(np.percentile(r, 95), decimals=1) for r in results_df['stats list']]
    results_df['5% p-val'] = [np.around(np.percentile(r, 5), decimals=2) for r in results_df['p-val list']]
    results_df['Median p-val'] = [np.around(r[0], decimals=2) for r in results_df['p-val list']]
    results_df['95% p-val'] = [np.around(np.percentile(r, 95), decimals=2) for r in results_df['p-val list']]
    # process information from score differences
    results_df['All mutation'] = [np.sum(np.array(r)!=0) for r in results_df['score difference']]
    results_df['Pos mutation'] = [np.sum(np.array(r)>0) for r in results_df['score difference']]
    results_df['Neg mutation'] = [np.sum(np.array(r)<0) for r in results_df['score difference']]
    median_diff_list = []
    mean_diff_list = []
    for r in results_df['score difference']:
        nonzero_diff = np.array(r)[np.array(r)!=0]
        if len(nonzero_diff) > 0:
            median_diff_list.append(np.median(nonzero_diff))
            mean_diff_list.append(np.mean(nonzero_diff))
        else:
            median_diff_list.append(0)
            mean_diff_list.append(0)
    results_df['Median score difference'] = np.around(median_diff_list, decimals=3)
    results_df['Mean score difference'] = np.around(mean_diff_list, decimals=3)
    # make the final output format
    results_df = results_df[['motif', '5% stats', 'Median stats', '95% stats', 
                             '5% p-val', 'Median p-val', '95% p-val', 
                             'All mutation', 'Pos mutation', 'Neg mutation',
                             'Median score difference', 'Mean score difference', 'score difference']]
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
        merge_contain = np.array(list(s))
        mean_set = input_df.loc[merge_contain, ['5% stats', 'Median stats', '95% stats', 
                                                '5% p-val', 'Median p-val', '95% p-val']].mean(axis=0)
        mean_set = np.around(mean_set, decimals=2)
        # sort merged motifs by their absolute log p-values
        sorted_idx = np.argsort(np.abs(np.array(input_df.loc[merge_contain, 'Median p-val'])))[::-1]
        top_summary = input_df.loc[merge_contain[sorted_idx[0]], 
                                   ['All mutation', 'Pos mutation', 'Neg mutation', 
                                    'Median score difference', 'Mean score difference']
                                  ]
        merge_contain = merge_contain[sorted_idx]
        merge_name = '|'.join([input_df.loc[ss, 'motif']+'$'+ss for ss in merge_contain])
        mean_set.name = merge_name
        top_summary.name = merge_name
        merge_stats = pd.concat([merge_stats, pd.concat([mean_set, top_summary], axis=0)], axis=1)
    merge_stats = merge_stats.T
    merge_stats.index.name = 'Merged motif (Total sequences: '+str(len(input_df['score difference'][0]))+')'
    merge_stats = merge_stats[['5% stats', 'Median stats', '95% stats', 
                               '5% p-val', 'Median p-val', '95% p-val', 
                               'All mutation', 'Pos mutation', 'Neg mutation', 
                               'Median score difference', 'Mean score difference']]

    return merge_stats    
    
