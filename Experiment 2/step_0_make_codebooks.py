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


if __name__ == '__main__':
    codebooks = []
    codons = list(codon_freq_dict.keys())

    # # Make codebooks with varying popularity
    p = np.array(list(codon_freq_dict.values()))
    for k in np.arange(-5,5.5,0.5):
        p_modified = p**k   # raising probability vector to higher powers makes popular codons even more popular.
        p_modified = p_modified/p_modified.sum()
        codebook = np.sort(np.random.choice(codons, 16, replace=False, p=p_modified))
        codebooks.append(codebook)
    for k in np.arange(-2,2,0.2):
        p_modified = p**k   # raising probability vector to higher powers makes popular codons even more popular.
        p_modified = p_modified/p_modified.sum()
        codebook = np.sort(np.random.choice(codons, 16, replace=False, p=p_modified))
        codebooks.append(codebook)
    pop_codeboooks = len(codebooks)
    print(f'Made {pop_codeboooks} codebooks for testing popularity')

    # Make codebooks with varying GC content
    gc_codons = np.zeros(len(codons))
    for i, codon in enumerate(codons):
        gc_codons[i] = codon.count('G') + codon.count('C') + 1      # gives a number from 1-4. 

    for k in np.arange(7.5,2.5,-0.5):
        p_low_gc = (5-gc_codons)**k
        p_low_gc = p_low_gc / p_low_gc.sum()
        codebook = np.sort(np.random.choice(codons, 16, replace=False, p=p_low_gc))
        codebooks.append(codebook)
    for k in np.arange(2.5,0,-0.2):
        p_low_gc = (5-gc_codons)**k
        p_low_gc = p_low_gc / p_low_gc.sum()
        codebook = np.sort(np.random.choice(codons, 16, replace=False, p=p_low_gc))
        codebooks.append(codebook)
    for k in np.arange(0,2.5,0.2):
        p_high_gc = gc_codons**k
        p_high_gc = p_high_gc / p_high_gc.sum()
        codebook = np.sort(np.random.choice(codons, 16, replace=False, p=p_high_gc))
        codebooks.append(codebook)
    for k in np.arange(2.5,7.5,0.5):
        p_high_gc = gc_codons**k
        p_high_gc = p_high_gc / p_high_gc.sum()
        codebook = np.sort(np.random.choice(codons, 16, replace=False, p=p_high_gc))
        codebooks.append(codebook)
    GC_codebooks = len(codebooks) - pop_codeboooks
    print(f'Made {GC_codebooks} codebooks for testing GC content')
    
    
    # Make codebooks with varying homopolymer metric
    # Cannot code for pos_31 2mers, so we focus on pos_12 and pos_23 2mers
    homopoly_codons = np.zeros(len(codons))
    for i, codon in enumerate(codons):
        homopoly_codons[i] = int(codon[0]==codon[1]) + int(codon[1]==codon[2]) + 1          # give a number from 1 to 3
    for k in np.arange(7.5,2.5,-0.5):
        p_low_2mer = (4-homopoly_codons)**k
        p_low_2mer = p_low_2mer / p_low_2mer.sum()
        codebook = np.sort(np.random.choice(codons, 16, replace=False, p=p_low_2mer))
        codebooks.append(codebook)
    for k in np.arange(2.5,0,-0.2):
        p_low_2mer = (4-homopoly_codons)**k
        p_low_2mer = p_low_2mer / p_low_2mer.sum()
        codebook = np.sort(np.random.choice(codons, 16, replace=False, p=p_low_2mer))
        codebooks.append(codebook)
    for k in np.arange(0,2.5,0.2):
        p_high_2mer = homopoly_codons**k
        p_high_2mer = p_high_2mer / p_high_2mer.sum()
        codebook = np.sort(np.random.choice(codons, 16, replace=False, p=p_high_2mer))
        codebooks.append(codebook)
    for k in np.arange(2.5,7.5,0.5):
        p_high_2mer = homopoly_codons**k
        p_high_2mer = p_high_2mer / p_high_2mer.sum()
        codebook = np.sort(np.random.choice(codons, 16, replace=False, p=p_high_2mer))
        codebooks.append(codebook)
    homo_codebooks = len(codebooks) - pop_codeboooks - GC_codebooks
    print(f'Made {homo_codebooks} codebooks for testing homopolymerity')

    # Add codebooks to folder
    for i, codebook in enumerate(codebooks):
        with open(f'./1_codebooks/codebook_{i+1}.txt', 'w') as f:
            codebook_string = ' '.join(codebook)
            f.write(codebook_string)