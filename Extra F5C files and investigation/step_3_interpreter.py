'''
This code will looks at the events tsv files and estimate a consensus squiggle via averaging. 

'''
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import csv
from itertools import product
from step_0_configure import read_length, kmer_len
import os
import sys



os.chdir(sys.path[0])

# Create dictionary with {kmer: (kmer_mean, kmer_std)} entries for all 4^9 possible 9mers.
# R10.model is copied from one of the models in the model.h file in the f5c source files. 
# This is distinct from the file found on kmer_models git page.
model10path_f5c = './R10.model'
model10f5c = {}



with open('consensus_squiggles_scaled.csv', 'w', newline='') as csv_output:
    writer = csv.writer(csv_output)
    for (direction, n) in product(['reverse', 'template'], range(1, 9)):
        sequence = f'{direction}_{n}'
        print(f'Analysing {sequence}')
        events = pd.read_csv(f'./5_aligned_events/events_{sequence}_scaled_id_collapse.tsv', delimiter='\t')

        # Mapping between f5c position columns and reference_kmer column
        position_to_kmer = events[['position', 'reference_kmer']].drop_duplicates().set_index('position')

        current_levels_weighted = [sequence]
        

        for i in range(read_length-kmer_len+1):
            # Look only at rows corresponding to current position/kmer. Remove rows with infinite standardised level.
            events_filt = events[events.position == i][['position','read_name','event_level_mean','event_stdv','event_length']].replace([np.inf, -np.inf], np.nan).dropna(axis=0, how='any')
            # Take duration/std-weighted mean of event current levels
            mu_vec = events_filt['event_level_mean']
            L_vec = events_filt['event_length']
            std_vec = np.maximum(events_filt['event_stdv'],1)
            current_levels_weighted.append( (mu_vec*L_vec/std_vec).sum() / (L_vec/std_vec).sum() )

        writer.writerow(current_levels_weighted)



# Older consensus method
'''
with open(model10path_f5c, 'r') as tsvfile:
    reader = csv.reader(tsvfile, delimiter=',')
    for row in reader:
        row = ([thing.replace('/','').strip() for thing in row])
        model10f5c[row[2]] = (float(row[0]), float(row[1]))

        
with open('consensus_squiggles.csv', 'w', newline='') as csv_output:
    writer = csv.writer(csv_output)
    for (direction, n) in product(['reverse', 'template'], range(1, 9)):
        sequence = f'{direction}_{n}'
        print(f'Analysing {sequence}')
        events = pd.read_csv(f'./5_aligned_events/events_{sequence}.tsv', delimiter='\t')

        # Mapping between f5c position columns and reference_kmer column
        position_to_kmer = events[['position', 'reference_kmer']].drop_duplicates().set_index('position')

        standardised_levels_duration_weighted = [sequence, 'duration weighted']
        current_levels_duration_weighted = [sequence, 'duration weighted']

        for i in range(read_length-kmer_len+1):
            # Look only at rows corresponding to current position/kmer. Remove rows with infinite standardised level.
            events_filt = events[events.position == i][['position','read_index','standardized_level','event_length']].replace([np.inf, -np.inf], np.nan).dropna(axis=0, how='any')
            # Take duration-weighted mean of standardized levels
            standardised_levels_duration_weighted.append((events_filt['standardized_level'] * events_filt['event_length']).sum() / (events_filt['event_length']).sum())
            reference_kmer = position_to_kmer.loc[i, 'reference_kmer']
            # Current average standardized level to a current level.
            current_levels_duration_weighted.append(float(model10f5c[reference_kmer][0]) + standardised_levels_duration_weighted[-1] * float(model10f5c[reference_kmer][1]))

        standardised_levels_read_weighted = [sequence, 'unweighted']
        current_levels_read_weighted = [sequence, 'unweighted']

        for i in range(read_length-kmer_len+1):
            # Look only at rows corresponding to current position/kmer. Remove rows with infinite standardised level.
            events_filt = events[events.position == i][['position','read_index','standardized_level','event_length']].replace([np.inf, -np.inf], np.nan).dropna(axis=0, how='any')
            standardised_levels_read_weighted.append(events_filt.groupby('read_index').apply(lambda x: (x['standardized_level'] * x['event_length']).sum() / x['event_length'].sum(), include_groups=False).mean())
            reference_kmer = position_to_kmer.loc[i, 'reference_kmer']
            current_levels_read_weighted.append(float(model10f5c[reference_kmer][0]) + standardised_levels_read_weighted[-1] * float(model10f5c[reference_kmer][1]))
        
        writer.writerow(standardised_levels_duration_weighted)
        writer.writerow(current_levels_duration_weighted)
        writer.writerow(standardised_levels_read_weighted)
        writer.writerow(current_levels_read_weighted)
'''