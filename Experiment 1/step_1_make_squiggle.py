"""
Experiment 1 Step 1
step_1_make_squiggles.py
First adjust parameters to desired values in step_0_configure.py, then run file
For each combination of parameters, a file will appear in 2_tx_msg/ and 3_squiggles/ 
"""
import random
import numpy as np
import time
import csv
from itertools import product
from step_0_configure import len_msg, PAD_str, num_trials, sd_min, sd_max, dwell_min, dwell_max, choice_list, mu, sigma, offset, scale, dorado_call_auto, kmer_len, kmer_centre, squiggle_generation, num_codebooks

import os
import sys
os.chdir(sys.path[0])

# Scrappie only works on Linux, so import is conditional
if squiggle_generation == 'scrappie':
    import scrappy
else:
    from step_0_generate_kmer_dict import kmer_dict

# Retrieve the four codebooks used in Experiment 1
codebooks = [None]*num_codebooks
for i in range(num_codebooks):
    with open(f'./1_codebooks/codebook_{i+1}.txt','r') as f:
        codebooks[i] = f.readlines()[0].split(' ')

# Labels describing each codebook
option_labels = [   'Most popular codons', 
                    'Least popular codons', 
                    'Random codons', 
                    'Hard-coded bad codons'] + list(range(5, num_codebooks+1))

# Function to save squiggles to CSV
def save_list_of_lists_to_csv(filename, data):
    with open(filename, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerows(data)


# Iterate over each combinations of trial parameters
total_time = 0
total_combinations = len(choice_list) * (sd_max+1-sd_min) * (dwell_max+1-dwell_min)
progress = 0
global_start_time = time.time()
for choice, sd_mean_scaled, dwell_mean in product(choice_list, range(sd_min,sd_max+1), range(dwell_min,dwell_max+1)):
    # Time the trials to estimate total runtime
    trial_start_time = time.time()

    # File names
    output_file = f'./2_tx_msg/tx_message_sd{sd_mean_scaled}_K{dwell_mean}_cb{choice}.fa'
    csv_filename = f'./3_squiggles/squiggles_sd{sd_mean_scaled}_K{dwell_mean}_cb{choice}.csv'

    # Obtain codebook and scale sd
    codebook = codebooks[choice-1]
    sd_mean = sd_mean_scaled / 100

    # Initialise trial
    description = f"Option {choice}: {option_labels[choice-1]}"
    print("#"*90)
    print(f"Description: {description}.\nConducting {num_trials} trials with sd={sd_mean:.2f} and dwell={dwell_mean}")
    all_signals = []
    transmit_text = ""

    # Conduct repeated trials
    for i in range(num_trials):
        # Generate transmission message
        tx_msg = ('').join(random.choices(codebook, k=len_msg))

        # Pad message with PAD_str
        tx_msg_padded = PAD_str + tx_msg + PAD_str

        # Create read id, which is @<sd_mean>-<dwell_mean>-<codebook>-0000-<trial_number>
        read_id = f"@{str(sd_mean_scaled).zfill(8)}-{str(dwell_mean).zfill(4)}-{str(choice).zfill(4)}-0000-{str(i).zfill(12)}"
        transmit_text += (read_id + '\n' + tx_msg + '\n')

        # Obtain squiggle statistics
        padded_len = len(tx_msg_padded)
        # Obtain squiggle statistics using scrappie
        if squiggle_generation == 'scrappie':
            scrappy_stats = scrappy.sequence_to_squiggle(tx_msg_padded, rescale=True)
            tx_msg_mean_values = scrappy_stats.data(as_numpy=True)[:,1]
        # Obtain squiggle statistics using kmer table
        elif squiggle_generation == 'kmer':
            tx_msg_mean_values = []
            # Head and tail are the number of bases before and after the centre of the kmer.
            head = kmer_centre - 1
            tail = kmer_len - kmer_centre
            for j in range(head, padded_len-tail):
                tx_msg_mean_values.append(kmer_dict[tx_msg_padded[j-head:j+tail+1]])
        
        # Generate sample values for each base
        signal_int16 = []
        for i in range(len(tx_msg_mean_values)):
            mean_curr = tx_msg_mean_values[i]
            # Dwell time follows a geometric distribution
            K_i = np.random.geometric(1/(dwell_mean)) 
            for _ in range(K_i):
                # Noise follows a Gaussian distribution
                current_norm = np.random.normal(loc=mean_curr, scale=sd_mean)
                # Normalised current level and pA current level can also be calculated.
                # signal_norm.append(signal_norm, current_norm)
                # signal_pa.append(signal_pa, sigma*current_norm+mu)
                signal_int16.append(np.uint16(round((sigma*current_norm+mu)/scale-offset)))
        all_signals.append(signal_int16)
    # Write all transmission messages to FASTA file
    with open(output_file, "w") as fasta_file:
        fasta_file.write(transmit_text)

    # Write all transmission squiggles to CSV file
    save_list_of_lists_to_csv(csv_filename, all_signals)

    # Print report of trial run time
    trial_time = time.time()-trial_start_time
    total_time += trial_time
    progress += 1
    est_time = total_time/progress * (total_combinations-progress)
    print(f'Trial set took {trial_time:.1f} seconds. Est. remaining time: {est_time:.1f} seconds')
    print("#"*90 + '\n')

print(f'Total squiggle generation time was {total_time:.1f} seconds')

if dorado_call_auto:
    import step_2_make_pod5_reads

print(f'Global run time was {time.time()-global_start_time:.2f} seconds')
