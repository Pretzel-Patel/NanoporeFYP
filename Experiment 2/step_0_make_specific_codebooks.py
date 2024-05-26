"""

"""

import pandas as pd
import numpy as np
import random
import os
import sys

os.chdir(sys.path[0])

# Codebook generation
# TSV file
tsv_file_location = './0_codon_statistics/nuclear_codon_statistics.tsv'

# Import popularity as dataframe
df = pd.read_csv(tsv_file_location, sep='\t')
df = df.sort_values(by='Frequency in 1000', ascending=False)
# Mapping codons to frequency
codon_freq_df = df[['CODON', 'Frequency in 1000']].set_index('CODON')
codon_freq_dict = (codon_freq_df/1000).squeeze().to_dict()


# Evaluate codebook and return three metrics
def evaluate_codebook(codebook):
    # Popularity metric - defined as average frequency of codebook divided by average frequency of all codons (1/64)
    popularity = 0
    for codon in codebook:
        popularity += codon_freq_dict[codon]
    popularity = popularity / len(codebook) * 64

    # GC content metric
    GC_content = 0
    for codon in codebook:
        GC_content += codon.count('G') + codon.count('C')
    GC_content = GC_content / (3*len(codebook))

    # Homopolymer metric
    # position 1&2
    pos_12 = 0
    pos_23 = 0
    pos_31 = 0
    for codon in codebook:
        if codon[0] == codon[1]:
            pos_12 += 1
        if codon[1] == codon[2]:
            pos_23 += 1
    for first_codon in codebook:
        for second_codon in codebook:
            if first_codon[2] == second_codon[0]:
                pos_31 += 1
    pos_12 = pos_12 / len(codebook)
    pos_23 = pos_23 / len(codebook)
    pos_31 = pos_31 / len(codebook)**2
    homopolymerity = (pos_12 + pos_23 + pos_31) / 3

    return (popularity, GC_content, homopolymerity)

def check_codebook_is_specific(codebook, metric):
    popularity, GC_content, homopolymerity = evaluate_codebook(codebook)
    if metric == 'popularity':
        return (0.46 < GC_content < 0.54) and (0.27 < homopolymerity < 0.33)
    elif metric == 'gc content':
        return (0.95 < popularity < 1.05) and (0.27 < homopolymerity < 0.33)
    elif metric == 'homopolymerity':
        return (0.95 < popularity < 1.05) and (0.46 < GC_content < 0.54)

if __name__ == '__main__':
    codebooks = []
    codons = list(codon_freq_dict.keys())


    # Make codebooks with varying popularity
    p = np.array(list(codon_freq_dict.values()))
    for k in np.arange(-5,5.5,0.1):
        while True:
            p_modified = p**k   # raising probability vector to higher powers makes popular codons even more popular.
            p_modified = p_modified/p_modified.sum()
            codebook = np.sort(np.random.choice(codons, 16, replace=False, p=p_modified))
            if check_codebook_is_specific(codebook, 'popularity'):
                break
        codebooks.append(codebook)
    for k in np.arange(-2,2,0.0425):
        while True:
            p_modified = p**k   # raising probability vector to higher powers makes popular codons even more popular.
            p_modified = p_modified/p_modified.sum()
            codebook = np.sort(np.random.choice(codons, 16, replace=False, p=p_modified))
            if check_codebook_is_specific(codebook, 'popularity'):
                break
        codebooks.append(codebook)

    pop_codeboooks = len(codebooks)
    print(f'Made {pop_codeboooks} codebooks for testing popularity')

    # Make codebooks with varying GC content
    gc_codons = np.zeros(len(codons))
    for i, codon in enumerate(codons):
        gc_codons[i] = codon.count('G') + codon.count('C') + 1      # gives a number from 1-4. 

    for k in np.arange(7.5,2.5,-0.1):
        while True:
            p_low_gc = (5-gc_codons)**k
            p_low_gc = p_low_gc / p_low_gc.sum()
            codebook = np.sort(np.random.choice(codons, 16, replace=False, p=p_low_gc))
            if check_codebook_is_specific(codebook, 'gc content'):
                break
        codebooks.append(codebook)
    for k in np.arange(2.5,0,-0.05):
        while True:            
            p_low_gc = (5-gc_codons)**k
            p_low_gc = p_low_gc / p_low_gc.sum()
            codebook = np.sort(np.random.choice(codons, 16, replace=False, p=p_low_gc))
            if check_codebook_is_specific(codebook, 'gc content'):
                break
        codebooks.append(codebook)
    for k in np.arange(0,2.5,0.05):
        while True:           
            p_high_gc = gc_codons**k
            p_high_gc = p_high_gc / p_high_gc.sum()
            codebook = np.sort(np.random.choice(codons, 16, replace=False, p=p_high_gc))
            if check_codebook_is_specific(codebook, 'gc content'):
                break
        codebooks.append(codebook)
    for k in np.arange(2.5,7.5,0.1):
        while True:           
            p_high_gc = gc_codons**k
            p_high_gc = p_high_gc / p_high_gc.sum()
            codebook = np.sort(np.random.choice(codons, 16, replace=False, p=p_high_gc))
            if check_codebook_is_specific(codebook, 'gc content'):
                break
        codebooks.append(codebook)
    GC_codebooks = len(codebooks) - pop_codeboooks
    print(f'Made {GC_codebooks} codebooks for testing GC content')
    
    
    # Make codebooks with varying homopolymer metric
    # Cannot code for pos_31 2mers, so we focus on pos_12 and pos_23 2mers
    homopoly_codons = np.zeros(len(codons))
    for i, codon in enumerate(codons):
        homopoly_codons[i] = int(codon[0]==codon[1]) + int(codon[1]==codon[2]) + 1          # give a number from 1 to 3
    for k in np.arange(7.5,2.5,-0.1):
        while True:
            p_low_2mer = (4-homopoly_codons)**k
            p_low_2mer = p_low_2mer / p_low_2mer.sum()
            codebook = np.sort(np.random.choice(codons, 16, replace=False, p=p_low_2mer))
            if check_codebook_is_specific(codebook, 'homopolymerity'):
                break
        codebooks.append(codebook)
    for k in np.arange(2.5,0,-0.05):
        while True:
            p_low_2mer = (4-homopoly_codons)**k
            p_low_2mer = p_low_2mer / p_low_2mer.sum()
            codebook = np.sort(np.random.choice(codons, 16, replace=False, p=p_low_2mer))
            if check_codebook_is_specific(codebook, 'homopolymerity'):
                break
        codebooks.append(codebook)
    for k in np.arange(0,2.5,0.05):
        while True:
            p_high_2mer = homopoly_codons**k
            p_high_2mer = p_high_2mer / p_high_2mer.sum()
            codebook = np.sort(np.random.choice(codons, 16, replace=False, p=p_high_2mer))
            if check_codebook_is_specific(codebook, 'homopolymerity'):
                break
        codebooks.append(codebook)
    for k in np.arange(2.5,7.5,0.1):
        while True:
            p_high_2mer = homopoly_codons**k
            p_high_2mer = p_high_2mer / p_high_2mer.sum()
            codebook = np.sort(np.random.choice(codons, 16, replace=False, p=p_high_2mer))
            if check_codebook_is_specific(codebook, 'homopolymerity'):
                break
        codebooks.append(codebook)
    homo_codebooks = len(codebooks) - pop_codeboooks - GC_codebooks
    print(f'Made {homo_codebooks} codebooks for testing homopolymerity')

    # Add codebooks to folder
    for i, codebook in enumerate(codebooks):
        with open(f'./1_codebooks/codebook_{i+1}.txt', 'w') as f:
            codebook_string = ' '.join(codebook)
            f.write(codebook_string)