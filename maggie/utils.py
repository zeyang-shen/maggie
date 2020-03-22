import numpy as np
import pandas as pd

import Bio
from Bio import motifs, SeqIO

def read_fasta(fasta_file, skip_duplicate=True):
    '''
    Read in sequences
    '''
    alphabet = Bio.Seq.IUPAC.Alphabet.IUPAC.IUPACUnambiguousDNA() # need to use this alphabet for motif score calculation
    id_seq_dict = {} # {sequenceID: fastq sequence}
    duplicate_keys = []
    for seq_record in SeqIO.parse(fasta_file, "fasta"):  
        seq_record.seq.alphabet = alphabet
        if seq_record.id in id_seq_dict.keys():
            duplicate_keys.append(seq_record.id)
        else:
            id_seq_dict[seq_record.id] = seq_record.seq
    # delete duplicate keys
    if skip_duplicate:
        for dk in duplicate_keys:
            del id_seq_dict[dk]
        if len(duplicate_keys) > 0:
            print('Ignore duplicate keys in %s: %s' % (fasta_file, duplicate_keys))
    return id_seq_dict


def FDR_cutoff(input_data, alpha=0.05):
    '''
    Find significant motif(s) based on corrected p-values after Benjamini–Hochberg procedure
    
    input_data:
        Require pandas series format with motif names as indices
    
    Parameters:
        alpha: cutoff for corrected p-values after Benjamini–Hochberg procedure (default: 0.05)
    
    Outputs:
        List of motifs passing significance threshold
    '''
    id_sorted = np.argsort(input_data)[::-1]
    pvalue_sorted = input_data[id_sorted]
    motif_sorted = input_data.index.values[id_sorted]
    m = len(pvalue_sorted)
    sort_i = 0
    for i in range(m):
        thr = -np.log10((i+1)*alpha/m)
        if pvalue_sorted[i] > thr:
            sort_i = i+1
    sig_motifs = motif_sorted[:sort_i]
    
    return sig_motifs


def write_fasta(seq_dict, out_path):
    '''
    Write a dictionary containing sequences to a fasta file
    '''
    string = ""
    for seq in seq_dict:
        string += ">" + str(seq) + "\n"
        string += str(seq_dict[seq]) + "\n"
    with open(out_path, 'w') as g:
        g.write(string)
    return


def cat_fasta_files(orig_fasta_files, mut_fasta_files):
    '''
    Concatenate sequences coming from multiple fasta files
    
    Input:
        orig_fasta_files, mut_fasta_files required to be aligned in pairs
    
    '''
    # combine sequences to dictionary
    orig_seq_dict = dict()
    mut_seq_dict = dict()
    for i in range(len(orig_fasta_files)):
        # original files
        tmp_dict1 = read_fasta(orig_fasta_files[i])
        tmp_dict1 = dict([(orig_fasta_files[i]+'|'+k, v) for k, v in tmp_dict1.items()])
        orig_seq_dict.update(tmp_dict1)
        # mutated files
        tmp_dict2 = read_fasta(mut_fasta_files[i])
        tmp_dict2 = dict([(mut_fasta_files[i]+'|'+k, v) for k, v in tmp_dict2.items()])
        mut_seq_dict.update(tmp_dict2)
    return orig_seq_dict, mut_seq_dict
    
### functions used in simulated data analysis
