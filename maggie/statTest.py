import numpy as np
import pandas as pd

import os

import Bio
from Bio import motifs
from Bio import SeqIO

from scipy.stats import ttest_1samp


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


def Ttest(X, motifs=None, save=True, prefix=None):
    X = np.array(X)
    p_values = np.array([ttest_1samp(i[i != 0], 0, axis=0)[1] for i in X])
    t_values = np.array([ttest_1samp(i[i != 0], 0, axis=0)[0] for i in X])
    p_values[np.isnan(p_values)] = 1
    t_values[np.isnan(t_values)] = 0
    
    if save:
        # filenames
        if prefix is None:
            save_file = './Ttest_stats.tsv'
        else:
            save_file = prefix + '_Ttest_stats.tsv'
        
        # save statistics
        if motifs is None:
            motifs = ['']*len(X)
        with open(save_file, 'w') as file:
            file.write('\t'.join(['motif', 'p-value', 't-value']) + '\n')
            for i in range(len(X)):
                file.write('\t'.join([motifs[i], str(p_values[i]), str(t_values[i])]) + '\n')  
    else:
        return p_values, t_values


def Filter(stats_files, sig_t=0.01, tpm_file=None, exp_t=5):
    # read statistics files
    all_pv_df = pd.DataFrame()
    for file in stats_files:
        stats_df = pd.read_csv(file, sep='\t', index_col=0)
        pv = stats_df['p-value']
        tv = stats_df['t-value']
        bool_pv = -np.log10(pv)*np.sign(tv)
        all_pv_df = pd.concat([all_pv_df, bool_pv], axis=1)
    all_pv_df.columns = stats_files
    
    # 1st filter by p-values
    p_thresh = -np.log10(sig_t/len(all_pv_df))
    pv_log_filt_df = all_pv_df.loc[np.array(np.abs(all_pv_df).max(axis=1) > p_thresh)]
    pv_log_filt_df = pv_log_filt_df.sort_values(stats_files[0], ascending=False)
    
    if tpm_file is None:
        return pv_log_filt_df, None
    
    # 2nd filter by expression level
    else:
        # read TPM file
        tpm_file_df = pd.read_csv(tpm_file, sep='\t', index_col=0)
        tpm_file_df.index = [i.split('|')[0] for i in tpm_file_df.loc[:,'Annotation/Divergence']]
        tpm_df = tpm_file_df.iloc[:,7:]
        tpm_log2_df = np.log2(tpm_df + 1)
        tpm_log2_df.index = [i.upper() for i in tpm_log2_df.index.values]
        
        # get TPM of significant motifs
        sig_list = [i.split('(')[0].upper().split('::') for i in pv_log_filt_df.index.values]
        sig_tpm = pd.DataFrame()
        no_expr = []
        for l in sig_list:
            try:
                sig_tpm = pd.concat([sig_tpm, tpm_log2_df.loc[l].min(axis=0)], axis=1)
            except:
                no_expr.extend(l)
                sig_tpm = pd.concat([sig_tpm, pd.Series([np.nan]*len(sig_tpm), index=sig_tpm.index)], axis=1)
        if len(no_expr) > 0:
            print('Missing gene expression:', no_expr)
        sig_tpm = sig_tpm.fillna(0)
        sig_tpm = sig_tpm.T
        sig_tpm.index = pv_log_filt_df.index
        
        # filter out low expressed motifs
        expr_bools = sig_tpm.max(axis=1)>exp_t
        expr_filt_df = pv_log_filt_df.loc[expr_bools]
        tpm_filt_df = sig_tpm.loc[expr_bools]
        return expr_filt_df, tpm_filt_df
    
    