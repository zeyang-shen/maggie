#!/usr/bin/env python
import os
import argparse

import sys
sys.path.append(os.path.dirname(__file__)+'/..')

from maggie import score, utils, visual, simulate

import numpy as np
import pandas as pd

def check_positive(value):
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue

def check_0_1_range(value):
    fvalue = float(value)
    if fvalue < 0 or fvalue > 1:
        raise argparse.ArgumentTypeError("%s is an invalid value" % value)
    return fvalue

if __name__ == "__main__":    
    parser = argparse.ArgumentParser(description='MAGGIE framework working with simulated datasets')
    parser.add_argument("insertMotif", 
                        help="motifs to insert into simulated sequences; multiple motifs should be separated by comma without space like 'SPI1,CEBPB'",
                        type=str)
    parser.add_argument("mutateMotifIndex", 
                        help="specify which motifs to be mutated; multiple indices should be separated by space like '1 2'. Index starts from 1",
                        nargs='+', type=int)
    parser.add_argument("sequenceNumber",
                        help="number of sequences to generate",
                        type=check_positive)
    parser.add_argument("-l", "--length",
                        help="length of sequence",
                        default=100,
                        type=check_positive)
    parser.add_argument("--saveSeq",
                        help="Flag for saving simulated sequences, default = False. Will generate two files that correspond to positive and negative sequences",
                        action='store_true')
    parser.add_argument("--mutation",
                        help="number of mutation/variant in each sequence",
                        default=1,
                        type=check_positive)
    parser.add_argument("--snr",
                        help="signal-to-noise ratio; ratio between target mutations and random mutations",
                        default=0.5,
                        type=check_0_1_range)
    parser.add_argument("--motifPath",
                        help="path to the motif files",
                        default=os.path.dirname(__file__)+"/../data/JASPAR2020_CORE_vertebrates_motifs/",
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
    
    insert_motifs = args.insertMotif
    mutation_index = args.mutateMotifIndex
    number = args.sequenceNumber
    length = args.length
    saveSeq = args.saveSeq
    snr = args.snr
    mutation = args.mutation
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
    
    # Read in motif files
    motif_dict = score.load_motifs(motif_dir)
    insert_motif_list = []
    for i in insert_motifs.split(','):
        found = False
        for k in motif_dict.keys():
            if k.upper() == i.upper() or k.split('$')[0].upper() == i.upper() or k.split('$')[1].upper() == i.upper():
                insert_motif_list.append(k)
                found = True
                break
        if not found:
            sys.exit('No matching motif found for '+str(i)+'. Please make sure all motifs are valid!')
    print('Inserted motifs:', insert_motif_list)
    new_mutation_index = []
    for mi in mutation_index:
        if mi <= 0 or mi > len(insert_motif_list):
            sys.exit('Invalid index for mutated motif: '+str(mi))
        else:
            new_mutation_index.append(mi-1)
    print('Mutated motifs:', np.array(insert_motif_list)[new_mutation_index])
        
    # Create folder to store outputs
    try:
        os.mkdir(output_dir)
    except FileExistsError:
        sys.exit('ERROR: Specified folder exists! Please name a different folder to save results')
    except OSError:
        sys.exit('ERROR: Please check if your specified path exists!')
    
    # Generate simulated sequences
    seq_list = simulate.simulate_sequence_pairs(length, number, 
                                                motif_dict, insert_motif_list, new_mutation_index, 
                                                mutation=mutation, snr=snr)
    pos_seqs = []
    neg_seqs = []
    for seq_tuple in seq_list:
        pos_seqs.append(seq_tuple[0])
        neg_seqs.append(seq_tuple[1])
    seq_df = pd.DataFrame(columns=['ref_seq', 'alt_seq'])
    seq_df['ref_seq'] = pos_seqs
    seq_df['alt_seq'] = neg_seqs
    orig_seq_dict = dict(seq_df['ref_seq'])
    mut_seq_dict = dict(seq_df['alt_seq'])
    if saveSeq:
        utils.write_fasta(orig_seq_dict, output_dir+'/positiveSeqs.fa')
        utils.write_fasta(mut_seq_dict, output_dir+'/negativeSeqs.fa')
        print('Successfully saved sequences of testing')
    
    # Run Maggie pipeline
    motif_list = list((motif_dict.keys()))
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
    