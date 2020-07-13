#!/usr/bin/env python
import sys
import os
import argparse
import requests

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download genomic data')
    parser.add_argument("genome", 
                        help="genome to download: hg19, hg38, mm10, hg18",
                        type=str)
    parser.add_argument("--annot", 
                        help="Flag for downloading annotation files at the same time, default = False",
                        action='store_true')
    parser.add_argument("-o", "--output", 
                        help="output directory; by default, will put the downloaded files in ./data/genomes/",
                        default=os.path.dirname(__file__)+'/../data/genomes/',
                        type=str)
    args = parser.parse_args()
    
    genome = args.genome
    annot = args.annot
    output_dir = args.output
    
    genome_path = output_dir+genome+'.fa'
    if os.path.exists(genome_path):
        valid = {"yes": True, "y": True, "ye": True, "no": False, "n": False}
        prompt = " [y/n] "
        while True:
            sys.stdout.write('WARNING: '+genome_path+' already exists! Are you sure to replace the current file?' + prompt)
            choice = input().lower()
            if choice in valid:
                if valid[choice]:
                    break
                else:
                    sys.exit('Terminated')
            else:
                sys.stdout.write("Please respond with 'yes' (or 'y') or 'no' (or 'n').\n")
            
    print('Downloading reference genome', genome)
    try:
        url = 'http://homer.ucsd.edu/zeyang/maggie/genomes/'+genome+'.fa'
        r = requests.get(url, stream=True)
        r.raise_for_status()
        total_size = int(r.headers.get('content-length', 0))
        block_size = 1024
        try:
            from tqdm import tqdm # load package to show the progress bar
            t=tqdm(total=total_size, unit='B', unit_scale=True)
            with open(genome_path, 'wb') as f:
                for data in r.iter_content(block_size):
                        t.update(len(data))
                        f.write(data)
                t.close()
        except:
            with open(genome_path, 'wb') as f:
                f.write(r.content)

    except:
        print('ERROR: failed to download the specified genome')
        raise
