## modules to generate simulated data
import numpy as np
import pandas as pd

from Bio import Seq


def generate_motif(motif_name, motif_dict):
    '''
    Randomly generate DNA sequence matching the PWM of any given motif
    '''
    nucleotide = ['A', 'C', 'G', 'T']
    motif_object = motif_dict[motif_name]
    length = len(motif_object)
    pwm = motif_object.pwm
    gen_motif = ''
    for i in range(length):
        pos = [0]
        for N in nucleotide:
            pos.append(pwm[N][i]+pos[-1])
        pos[-1] = 1
        n = np.random.random()
        ind = np.where(np.array(pos) < n)[0][-1]
        gen_motif += nucleotide[ind]
        
    return gen_motif


def simulate_sequence_pairs(length, number, motif_dict, insert_motif_list, mutation_index, mutation=1, snr=1):
    '''
    Generate simulated DNA sequence pairs with motif mutations
    
    Parameters:
        length: length of DNA sequence
        number: number of simulated sequences
        motif_dict: dictionary storing Biopython motifs using "load_motifs" module in "score.py"
        insert_motif_list: keys of motifs to be inserted in simulated sequences
        mutation_index: indexes of motifs in insert_motif_list to be mutated
        mutation: number of positions to be mutated
        snr: signal-noise-ratio, or proportion of simualted sequences experiencing mutations
    
    Outputs:
        A list of (positive sequence, negative sequence)
        
    Example:
        motif_dict = score.load_motifs('./data/JASPAR2020_CORE_vertebrates_motifs/')
        # Insert IRF1 and ZNF410 motifs and mutate only ZNF410 motif by 1 nucleotide for half of the sequences
        seq_list = simulate_sequence_pairs(100, 20, motif_dict, ['IRF1$MA0050.1', 'ZNF410$MA0752.1'], 1, mutation=1, snr=0.5)
    '''
    
    if type(insert_motif_list) is not list:
        insert_motif_list = [insert_motif_list]
    if type(mutation_index) is not list:
        mutation_index = [mutation_index]

    cut = length//len(insert_motif_list)
    if cut < 20:
        sys.exit('Too many motifs and not enough length of sequence')
    alphabet = Seq.IUPAC.Alphabet.IUPAC.IUPACUnambiguousDNA()
    seq_list = []
    for sn in range(number):
        seq = np.random.choice(['A', 'C', 'G', 'T'], replace=True, size=length)
        cut = length//len(insert_motif_list)
        insert_pos_list = []
        for o, im in enumerate(insert_motif_list):
            # generate a motif sequence
            motif_seq = generate_motif(im, motif_dict)
            # generate insertion position
            insert_pos = np.random.randint(cut*o, cut*(o+1)-len(motif_seq))
            insert_pos_list.append(insert_pos)
            # insert motif
            seq[insert_pos:insert_pos+len(motif_seq)] = list(motif_seq)

        # generate motif mutation
        mut_seq = seq.copy()
        for o in mutation_index:
            if sn < number*snr: #replace nucleotides of current motif sequence
                motif_len = len(motif_dict[insert_motif_list[o]])
                mutate_pos = np.random.choice(np.arange(insert_pos_list[o], insert_pos_list[o]+motif_len), 
                                              replace=False, size=mutation)
                for k in mutate_pos:
                    mut_seq[k] = np.random.choice(list(set(['A', 'C', 'G', 'T'])-set([mut_seq[k]])))
            else: #replace with another random motif sequence
                motif_seq = generate_motif(insert_motif_list[o], motif_dict)
                mut_seq[insert_pos_list[o]:insert_pos_list[o]+len(motif_seq)] = list(motif_seq)
        seq_list.append((''.join(seq), ''.join(mut_seq)))
        
    return seq_list