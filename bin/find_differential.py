#!/usr/bin/env python
import sys
import os
import argparse

import pandas as pd

def find_overlap(start1, end1, start2, end2):
    "find the overlap between the range (start1, end1) and the range (start2, end2)"
    return max(max((end2-start1), 0) - max((end2-end1), 0) - max((start2-start1), 0), 0)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find differential ranges from two BED files')
    parser.add_argument("-i1", 
                        help="BED file 1",
                        type=str)
    parser.add_argument("-i2", 
                        help="BED file 2",
                        type=str)
    parser.add_argument("-o1", 
                        help="output file to store sites specific to BED file 1",
                        type=str)
    parser.add_argument("-o2", 
                        help="output file to store sites specific to BED file 2",
                        type=str)
    parser.add_argument("--format", 
                        help="specify the output file format; default is BED",
                        choices=['peak', 'bed'])
    
    args = parser.parse_args()
    
    input_file1 = args.i1
    input_file2 = args.i2
    output_file1 = args.o1
    output_file2 = args.o2
    output_format = args.format
    
    #read in data
    d1 = pd.read_csv(input_file1, sep='\t', index_col=3, header=None)
    d2 = pd.read_csv(input_file2, sep='\t', index_col=3, header=None)
    chr_keys = set(d1[0])
    chr_keys.update(set(d2[0]))

    #find differential sites
    overlap_cutoff = 0
    score_cutoff = 0
    d1_spec = []
    d2_spec = []
    common = []
    for k in chr_keys:
        d1_cur = d1.loc[d1[0] == k]
        d2_cur = d2.loc[d2[0] == k]
        d1_cur = d1_cur.sort_values(by=1, ascending=True)
        d2_cur = d2_cur.sort_values(by=1, ascending=True)
        d1_idx = 0
        d2_idx = 0
        while d1_idx < len(d1_cur) and d2_idx < len(d2_cur):
            start1 = d1_cur.iloc[d1_idx, 1]
            end1 = d1_cur.iloc[d1_idx, 2]
            start2 = d2_cur.iloc[d2_idx, 1]
            end2 = d2_cur.iloc[d2_idx, 2]
            if start2 > end1:
                if d1_cur.iloc[d1_idx, 3] >= score_cutoff:
                    d1_spec.append(d1_cur.index[d1_idx])
                d1_idx += 1
                continue
            if start1 > end2:
                if d2_cur.iloc[d2_idx, 3] >= score_cutoff:
                    d2_spec.append(d2_cur.index[d2_idx])
                d2_idx += 1
                continue
            ov = find_overlap(start1, end1, start2, end2)
            if ov > overlap_cutoff:
                common.append((d1_cur.index[d1_idx], d2_cur.index[d2_idx], ov))
                d1_idx += 1
                d2_idx += 1
        if d1_idx < len(d1_cur):
            for rest_ in range(d1_idx, len(d1_cur)):
                d1_spec.append(d1_cur.index[rest_])

        if d2_idx < len(d2_cur):
            for rest_ in range(d2_idx, len(d2_cur)):
                d2_spec.append(d2_cur.index[rest_])
        
    #output differential sites
    d1_spec_df = d1.loc[d1_spec]
    d2_spec_df = d2.loc[d2_spec]
    if output_format == 'peak':
        d1_spec_df = d1_spec_df[[0,1,2,5]]
        d2_spec_df = d2_spec_df[[0,1,2,5]]
        d1_spec_df.index.name = 'PeakID'
        d2_spec_df.index.name = 'PeakID'
        d1_spec_df.columns = ['chr', 'start', 'end', 'strand']
        d2_spec_df.columns = ['chr', 'start', 'end', 'strand']
        d1_spec_df.to_csv(output_file1, sep='\t')
        d2_spec_df.to_csv(output_file2, sep='\t')
    elif output_format == 'bed':
        d1_spec_df[6] = d1_spec_df.index
        d2_spec_df[6] = d2_spec_df.index
        d1_spec_df.to_csv(output_file1, sep='\t', index=False, header=False)
        d2_spec_df.to_csv(output_file2, sep='\t', index=False, header=False)
    else:
        print('ERROR: Unexpected output format! Please specify either "peak" or "bed"')