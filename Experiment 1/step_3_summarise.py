import editdistance
import numpy as np
import math
import time
from itertools import product
import csv
import pyttsx3
import os
import sys
os.chdir(sys.path[0])

from step_0_configure import num_trials, sd_min, sd_max, dwell_min, dwell_max, choice_list, dorado_call_auto, model

def load_list_of_lists_from_csv(filename):
    with open(filename, 'r') as csvfile:
        csv_reader = csv.reader(csvfile)
        data = [list(map(int, row)) for row in csv_reader]
    return data



# Initialise dictionary containing (key, value) = (read_id, [tx_msg, rx_msg])
reads_dict = {}
for choice, sd_mean_scaled, dwell_mean in product(choice_list, range(sd_min,sd_max+1), range(dwell_min,dwell_max+1)):
    # File names
    tx_msg_file = f'./2_tx_msg/tx_message_sd{sd_mean_scaled}_K{dwell_mean}_cb{choice}.fa'
    # Extract text from each tx file
    with open(tx_msg_file,'r') as f:
        tx_msg_all = f.readlines()
    for i in range(num_trials):
        read_id = tx_msg_all[2*i][1:-1]         # strip @ and \n
        reads_dict[read_id] = [tx_msg_all[2*i+1][:-1],None]

# Either run basecalling with scrappie
total_time = 0
if model == 'scrappie':
    import scrappy
    for choice, sd_mean_scaled, dwell_mean in product(choice_list, range(sd_min,sd_max+1), range(dwell_min,dwell_max+1)):
        # File names
        print(f'Scrappie raw: codebook {choice}, sd {sd_mean_scaled/100}, K {dwell_mean}')
        csv_filename = f'./3_squiggles/squiggles_sd{sd_mean_scaled}_K{dwell_mean}_cb{choice}.csv'
        all_signals = load_list_of_lists_from_csv(csv_filename)
        for i in range(num_trials):
            signal_int16 = np.array(all_signals[i])
            read_id = f"0000{str(sd_mean_scaled).zfill(4)}-00{str(dwell_mean).zfill(2)}-00{str(choice).zfill(2)}-0000-0000{str(i).zfill(8)}"
            start = time.time()
            rx_msg = scrappy.basecall_raw(signal_int16)[0]          # Normalisation of raw_data is irrelevant.'
            total_time += time.time()-start
            reads_dict[read_id][1] = rx_msg
    print(f'Scrappy basecalling finished, taking {total_time:.2f} seconds.')
# or run basecalling with Dorado
else:
    if dorado_call_auto:
        with open(f'./5_rx_msg/all_calls.fa','r') as f:
            rx_msg_all = f.readlines()
    else:
        with open(f'./5_rx_msg/all_calls.fa','r',encoding='utf-16le') as f:       # Use this when command line dorado is used
            rx_msg_all = f.readlines()

    i = 0
    # rx_msg_all = [line.replace('\x00','') for line in rx_msg_all]
    valid_keys = []
    while i < len(rx_msg_all)-1:
        # print(i)
        read_id = (rx_msg_all[i].split('@'))[1][:-1]
        reads_dict[read_id][1] = rx_msg_all[i+1][:-1]
        i += 4


# Calculate edit distance between tx_msg and rx_msg
start = time.time()
results_dict = {}
edit_distances = []
for key in reads_dict:
    if reads_dict[key][1] == None:
        continue
    edit_dist = editdistance.eval(reads_dict[key][0], reads_dict[key][1])
    if key[0:18] in results_dict:
        results_dict[key[0:18]].append(edit_dist)
    else:
        results_dict[key[0:18]] = [edit_dist]
# Calculate length of transmission message
len_msg_bases = len(reads_dict[key][0])
print(f'Edit distance calculations took {time.time()-start:.2}s\n')

summary_all = []

for params, edit_distances in results_dict.items():
    sd_mean = int(params[0:8]) / 100
    dwell_mean = int(params[9:13])
    choice = int(params[14:])

    num_valid_trials = len(edit_distances)
    mean_result = np.mean(edit_distances)
    sd_result = np.std(edit_distances)
    Q_result = -10*math.log10(np.mean(edit_distances)/len_msg_bases)
    sd_pop = sd_result / np.sqrt(num_valid_trials)
    edit_upper = mean_result + 1.96*sd_pop
    edit_lower = mean_result - 1.96*sd_pop
    Q_upper = -10*math.log10(edit_lower/len_msg_bases)
    Q_lower = -10*math.log10(edit_upper/len_msg_bases)
    print(f"Parameters: codebook={choice}, sd={sd_mean}, K={dwell_mean}, num-trials={num_trials}, num-valid-trials={num_valid_trials}")
    print(f"Mean edit distance is {mean_result:.2f}, sd is {sd_result:.2f}, Q-score is {Q_result:.2f}")
    print(f"95% confidence interval for mean is [{edit_lower:.2f}, {edit_upper:.2f}] and for Q-score is [{Q_lower:.2f}, {Q_upper:.2f}]\n")
    excel_summary = f"{choice},{len_msg_bases//3},{num_valid_trials},{sd_mean},{dwell_mean},0,{mean_result:.2f},{sd_result:.2f}"
    # print(f"Excel summary: {excel_summary}")
    summary_all.append(excel_summary)


# Output summary of trials
output_file = "summary.txt"
with open(output_file, "w") as f:
    for line in summary_all:
        f.write(line)
        f.write('\n')

engine = pyttsx3.init()
engine.say('Trial is complete. Trial is complete. Trial is complete.')
engine.runAndWait()