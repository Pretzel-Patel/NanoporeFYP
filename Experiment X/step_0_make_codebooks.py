"""
Pipeline 4 is an adaptation of pipeline 2
Pipeline 4 runs multiple tests
"""

import pandas as pd
import numpy as np
import random
import os
import sys

os.chdir(sys.path[0])

# Process 1. Codebook generation
# TSV file
tsv_file_location = './0_codon_statistics/nuclear_codon_statistics.tsv'

# Option 1: Most common overall
df = pd.read_csv(tsv_file_location, sep='\t')
df = df.sort_values(by='Frequency in 1000', ascending=False)
# Choose 16 most frequent codons overall
option_1 = (list(df['CODON'][0:16]))
option_1.sort()

# Option 2: Least common overall
option_2 = (list(df['CODON'][-16:]))
option_2.sort()

# Option 3: Hard-coded fixed set
option_3 = ['AAT', 'AGT', 'ATA', 'ATT', 'CAA', 'CCC', 'CGG', 'CTC', 'CTT', 'GAG', 'GCC', 'GGA', 'GTG', 'TCG', 'TCT', 'TGT']

# Option 4: Hard-coded bad set
option_4 = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'AGA', 'ATA', 'CAA', 'GAA', 'TAA', 'TTT', 'TTA', 'TTC', 'TTG', 'TAT', 'TCT']

# Option 5: Random
bases = ['A', 'C', 'G', 'T']
option_5 = []
while len(option_5) < 16:
    candidate_3mer = ('').join(random.sample(bases,3))
    if (candidate_3mer in ['AAA', 'CCC', 'GGG', 'TTT']) or (candidate_3mer in option_5):
        continue
    option_5.append(candidate_3mer)
option_5.sort()

options = [option_1, option_2, option_3, option_4, option_5]
option_labels = [   'Most popular codons', 
                    'Least popular codons', 
                    'Random codons', 
                    'Hard-coded bad codons', 
                    'Random 3-mers', 
                    'External codebook']

for i in range(len(options)):
    print(option_labels[i])
    print(str(options[i]).replace("', '", " "))
    print()


'''
Archived
# Popular codons with no amino acid represnted twice
codon_freq_list = []
for amino_acid in set(df['Amino acid']):
    max_freq = 0
    max_codon = None
    for index, row in df.iterrows():
        if row['Amino acid'] == amino_acid and row['Frequency in 1000'] > max_freq:
            max_codon   = row["CODON"]
            max_freq    = row['Frequency in 1000']
    codon_freq_list.append((max_codon, max_freq))
# Remove homopolymers
for codon in codon_freq_list:
    if codon[0] in ['AAA', 'CCC', 'GGG', 'TTT']:
        codon_freq_list.remove(codon)
# Choose 16 most frequent codons, max 1 for each amino acid
while len(codon_freq_list) > 16:
    min_codon = (None,1000)
    for codon in codon_freq_list:
        if codon[1] < min_codon[1]:
            min_codon = codon
    codon_freq_list.remove(min_codon)
option_2 = []
for codon in codon_freq_list:
    option_2.append(codon[0])
option_2.sort()


'''