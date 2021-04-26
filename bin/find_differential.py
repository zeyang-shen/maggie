#!/usr/bin/env python
import sys
import os
import argparse
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find differential ranges from two BED files; bedtools required')
    parser.add_argument("-i1", 
                        help="BED file 1 (required)",
                        type=str)
    parser.add_argument("-i2", 
                        help="BED file 2 (required)",
                        type=str)
    parser.add_argument("-o1", 
                        help="output file to store sites specific to BED file 1 (optional)",
                        type=str)
    parser.add_argument("-o2", 
                        help="output file to store sites specific to BED file 2 (optional)",
                        type=str)
    parser.add_argument("--format", 
                        help="specify the output file format: HOMER peak file or bed file; default is BED",
                        choices=['peak', 'bed'], 
                        default='bed')
    
    args = parser.parse_args()
    
    input_file1 = args.i1
    input_file2 = args.i2
    output_file1 = args.o1
    output_file2 = args.o2
    output_format = args.format
    
    #set default output file
    if output_file1 is None:
        gettf1="echo `realpath "+str(input_file1)+ "`"
        gettf2="echo `realpath "+str(input_file2)+ "`"
        input_file1_path = os.popen(gettf1).read().strip()
        input_file2_path = os.popen(gettf2).read().strip()
        if output_format == 'bed':
            output_file1 = '/'.join(input_file1_path.split('/')[:-1])+'/input1_specific.bed'
            output_file2 = '/'.join(input_file1_path.split('/')[:-1])+'/input2_specific.bed'
        elif output_format == 'peak':
            output_file1 = '/'.join(input_file1_path.split('/')[:-1])+'/input1_specific.txt'
            output_file2 = '/'.join(input_file1_path.split('/')[:-1])+'/input2_specific.txt'
    
    #find differential sites
    cmd="bedtools intersect -a "+input_file1+" -b "+input_file2+" -v > "+output_file1
    result = os.system(cmd)
    if result != 0:
        sys.exit('Termination: failed to run bedtools!')
    cmd="bedtools intersect -a "+input_file2+" -b "+input_file1+" -v > "+output_file2
    result = os.system(cmd)
    if result != 0:
        sys.exit('Termination: failed to run bedtools!')
    
    #format to peak files
    if output_format == 'peak':
        d1 = pd.read_csv(output_file1, sep='\t', header=None)
        d2 = pd.read_csv(output_file2, sep='\t', header=None)
        d1 = d1[[3,0,1,2,5,4]]
        d2 = d2[[3,0,1,2,5,4]]
        d1.columns = ['ID', 'chr', 'start', 'end', 'strand', 'score']
        d2.columns = ['ID', 'chr', 'start', 'end', 'strand', 'score']
        d1.to_csv(output_file1, sep='\t', index=False, header=True)
        d2.to_csv(output_file1, sep='\t', index=False, header=True)
    print('Complete! Outputs saved in', output_file1, ',', output_file2)