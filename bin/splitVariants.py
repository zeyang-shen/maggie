#!/usr/bin/env python
import sys
sys.path.append('.')

from maggie import utils

import numpy as np
import pandas as pd

import os
import argparse
import requests

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Split variants based on genomic annotations')
    parser.add_argument("vcfFile", 
                        help="VCF file that contains testing variants",
                        type=str)
    parser.add_argument("genome", 
                        help="reference genome",
                        type=str)
    parser.add_argument("-L", "--overlap", 
                        help="overlap size to count for annotation",
                        type=int)
    parser.add_argument("-o", "--output", 
                        help="output directory; by default, will create a new folder under current path",
                        default='.',
                        type=str)
    parser.add_argument("-p", 
                        help="number of processors to run",
                        default=1,
                        type=int)
    args = parser.parse_args()
    
    vcf_file = args.vcfFile
    genome = args.genome
    overlap = args.overlap
    output_dir = args.output
    proc = args.p
    
    # Create folder to store outputs
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    
    # read in VCF file
    vcf_df = pd.read_csv(vcf_file, sep='\t', comment='#')
    if len(np.unique(vcf_df.iloc[:,2])) > len(vcf_df)//2:
        vcf_df.index = vcf_df.iloc[:,2]
    
    # download annotation file
    if len(genome.split('.')) == 1: # download reference genome
        gff_path = './data/genomes/'+genome+'.gff3'
        if not os.path.exists(gff_path):
            print('Downloading annotation file for', genome)
            try:
                url = 'http://homer.ucsd.edu/zeyang/maggie/genomes/'+genome+'.gff3'
                r = requests.get(url, stream=True)
                r.raise_for_status()
                total_size = int(r.headers.get('content-length', 0))
                block_size = 1024
                try:
                    from tqdm import tqdm # load package to show the progress bar
                    t=tqdm(total=total_size, unit='B', unit_scale=True)
                    with open(gff_path, 'wb') as f:
                        for data in r.iter_content(block_size):
                                t.update(len(data))
                                f.write(data)
                        t.close()
                except:
                    with open(gff_path, 'wb') as f:
                        f.write(r.content)
            except:
                print('ERROR: failed to download annotation file for the specified genome')
                raise
    else:
        gff_path = genome
    
    # Annotate regions
    size = 100
    region_df = pd.DataFrame(index=vcf_df.index, columns=['Chr', 'Start', 'End'])
    region_df['Chr'] = [str(i) for i in vcf_df.iloc[:,0]]
    region_df['Start'] = [int(i)-size//2 for i in vcf_df.iloc[:,1]]
    region_df['End'] = [int(i)+size//2 for i in vcf_df.iloc[:,1]]
    annotated_df = utils.annotateBED(gff_path, region_df, overlap=overlap)
    annotated_df.to_csv(output_dir+'/allAnnotation.txt', sep='\t')
    
    # Split regions based on annotations and save to separate files
    for annot_type in np.unique(annotated_df['Annotation']):
        split_vcf = vcf_df.loc[annotated_df.loc[annotated_df['Annotation'] == annot_type].index]
        split_vcf.to_csv(output_dir+'/variants_'+annot_type+'.vcf', sep='\t', index=None)
    