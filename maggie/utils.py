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


def load_genome(ref_path):
    '''
    Load reference genomes from a single fasta file
    '''
    ref_dict = {}
    for seq in SeqIO.parse(ref_path, "fasta"):
        chromID = seq.id
        chromSeq = (str(seq.seq)).upper()
        ref_dict[chromID] = chromSeq
    return ref_dict


def data_prep(path, genomes, size=100, skiprows=0, file_format='bed'):
    '''
    Extract sequences based on reference genome and input file with information of location
    '''
    data_list = []
    i = skiprows
    for line in open(path):
        # skip rows
        if i > 0:
            i -= 1
            continue
        
        # read each line
        elems = line.split()
        if file_format == 'bed':
            chromID = elems[0]
            start, end = int(elems[1]), int(elems[2])
        elif file_format == 'homer':
            chromID = elems[1]
            start, end = int(elems[2]), int(elems[3])
        elif file_format == 'vcf':
            chromID = elems[0]
            start, end = int(elems[1]), int(elems[1])
        
        # rescale the regions
        mid = (start+end)//2
        start = mid - size//2
        end = mid + size//2
        try:
            seq = genomes[chromID][start:end]
        except:
            seq = 'N'*size
        data_point = (seq, chromID, start, end)
        data_list.append(data_point)
    return data_list


def save_raw_data(results_df, output_path):
    '''
    Save score differences
    '''
    save_np = np.array([s for s in results_df['score difference']])
    np.save(output_path, save_np)
    print('Successfully saved score differences')


def annotateBED(gff_file, bed_df, overlap=None):
    '''
    Annotate regions in BED file based on GFF file
    '''
    # Read GFF file
    print('Reading annotation file')
    gencode = pd.read_table(gff_file, comment='#', sep='\t', 
                            names=['seqname', 'source', 'feature', 'start' , 'end', 
                                   'score', 'strand', 'frame', 'attribute'])
    
    combine_df = pd.DataFrame()
    for chrom in np.unique(bed_df.iloc[:,0]): # process chromosome-by-chromosome
        # obtain input regions
        chr_region = bed_df.loc[bed_df['Chr']==chrom].sort_values(by='Start')
        # obtain all annotations
        chr_gencode = gencode.loc[gencode['seqname']==chrom]
        if len(chr_gencode) == 0:
            continue
        # obtain transcripts
        chr_transcripts = chr_gencode.loc[chr_gencode['feature'] == 'transcript']
        chr_transcripts.is_copy = None
        chr_transcripts['ID'] = chr_transcripts.index.values
        chr_transcripts.index = [a.split('gene_name=')[1].split(';')[0] for a in chr_transcripts['attribute']]
        chr_transcripts = chr_transcripts.sort_values(by='start')
        chr_genes = chr_transcripts.loc[~chr_transcripts.index.duplicated(keep='first')]
        # filter genes
        gene_filter = []
        for g in chr_genes.index.values:
            a = chr_genes.loc[g, 'attribute']
            gene_type = a.split('gene_type=')[1].split(';')[0]
            if gene_type == 'protein_coding':
                gene_filter.append(g)
        chr_genes = chr_genes.loc[gene_filter]
        chr_transcripts = chr_transcripts.loc[gene_filter]
        print(chrom, len(chr_region), len(chr_gencode), len(chr_genes))
        
        # process region-by-region
        annots = []
        for row in chr_region.iterrows():
            bed_start, bed_end = row[1][1:3]
            snp_len = bed_end - bed_start
            snp_loc = (bed_start + bed_end)//2
            if overlap is None:
                overlap = snp_len//2

            # Get annotations for the nearest gene(s)
            near_ts_list = np.where(np.array(bed_end >= chr_transcripts['start']))[0]
            if len(near_ts_list) > 0:
                near_ts_idx = near_ts_list[-1]
                near_ts = chr_transcripts.iloc[near_ts_idx]
            else:
                near_ts_idx = 0
                near_ts = chr_transcripts.iloc[near_ts_idx]
            near_gene_idx = np.where(chr_genes.index==near_ts.name)[0][0]

            # Get annotations for the most confident transcript
            all_rela_ts = chr_transcripts.loc[[near_ts.name]]
            tsl = np.array([int(a.split('transcript_support_level=')[1].split(';')[0].replace('NA', '6')) 
                            if 'transcript_support_level' in a else 6 for a in all_rela_ts['attribute']])
            if min(tsl) <= 2:
                keep_ts_idx = np.where(tsl <= 2)[0]
            else:
                keep_ts_idx = np.argsort(tsl)[[0]]
            gene_comp = pd.DataFrame()
            for ti in keep_ts_idx:
                tmp_comp = chr_gencode.loc[all_rela_ts.iloc[ti]['ID']:]
                next_gene = np.where(np.any([tmp_comp['feature'][1:] == 'gene', 
                                             tmp_comp['feature'][1:] == 'transcript'], axis=0))[0]
                if len(next_gene) > 0:
                    tmp_comp = tmp_comp.iloc[:next_gene[0]]
                gene_comp = pd.concat([gene_comp, tmp_comp])

            # compute nearest gene and distance to TSS
            if near_gene_idx < len(chr_genes)-1:
                near_gene2 = chr_genes.index[near_gene_idx+1]
                near_gene2_all_ts = chr_transcripts.loc[[near_gene2]]
                tsl2 = np.array([int(a.split('transcript_support_level=')[1].split(';')[0].replace('NA', '6')) 
                                if 'transcript_support_level' in a else 6 for a in near_gene2_all_ts['attribute']])
                if min(tsl) <= 2:
                    keep_ts_idx2 = np.where(tsl2 <= 2)[0]
                else:
                    keep_ts_idx2 = np.argsort(tsl2)[[0]]
                near_gene2_ts = near_gene2_all_ts.iloc[keep_ts_idx2]
            all_near_ts = pd.concat([all_rela_ts.iloc[keep_ts_idx], near_gene2_all_ts.iloc[keep_ts_idx2]])
            all_dist_tss = []
            for ts in all_near_ts.iterrows():
                if ts[1]['strand'] == '+':
                    dist_tss = snp_loc - ts[1]['start']
                else:
                    dist_tss = snp_loc - ts[1]['end']
                all_dist_tss.append(dist_tss)
            nearest_idx = np.argsort(np.absolute(all_dist_tss))[0]
            dist_tss = all_dist_tss[nearest_idx]
            near_gene_name = all_near_ts.index[nearest_idx]

            # Assgin annotation
            # inside the transcript region
            if bed_start <= near_ts['end']-overlap:
                # 5 prime UTR
                gene_5UTR = gene_comp.loc[gene_comp['feature'] == 'five_prime_UTR']
                near_5UTR_idx = np.where(np.array(bed_end >= gene_5UTR['start']))[0]
                if len(near_5UTR_idx) > 0:
                    near_5UTR = gene_5UTR.iloc[near_5UTR_idx[-1]]
                    if bed_start <= near_5UTR['end']-overlap:
                        annots.append(('5\'UTR', dist_tss, near_ts.name, near_ts['strand']))
                        continue
                # 3 prime UTR
                gene_3UTR = gene_comp.loc[gene_comp['feature'] == 'three_prime_UTR']
                near_3UTR_idx = np.where(np.array(bed_end >= gene_3UTR['start']))[0]
                if len(near_3UTR_idx) > 0:
                    near_3UTR = gene_3UTR.iloc[near_3UTR_idx[-1]]
                    if bed_start <= near_3UTR['end']-overlap:
                        annots.append(('3\'UTR', dist_tss, near_ts.name, near_ts['strand']))
                        continue
                # coding region
                gene_exon = gene_comp.loc[gene_comp['feature'] == 'exon']
                near_exon_idx = np.where(np.array(bed_end >= gene_exon['start']))[0]
                if len(near_exon_idx) > 0:
                    near_exon = gene_exon.iloc[near_exon_idx[-1]]
                    if bed_start <= near_exon['end']-overlap:
                        annots.append(('exon', dist_tss, near_ts.name, near_ts['strand']))
                        continue   
                # if none of above, labeled as intron
                annots.append(('intron', dist_tss, near_ts.name, near_ts['strand']))

            # outside the transcript region but close to TSS
            elif np.abs(dist_tss) <= 1000:
                annots.append(('TSS', dist_tss, near_gene_name, near_ts['strand']))
            # outside the transcript region
            else:
                annots.append(('Intergenic', dist_tss, near_gene_name, near_ts['strand']))

        # Aggregate annotations
        annots = np.array(annots)
        annots_df = pd.DataFrame(annots, index=chr_region.index)
        annots_df.columns = ['Annotation', 'Distance to TSS', 'Nearest gene', 'Gene strand']
        combine_df = pd.concat([combine_df, pd.concat([chr_region, annots_df], axis=1)], axis=0)    
    return combine_df


