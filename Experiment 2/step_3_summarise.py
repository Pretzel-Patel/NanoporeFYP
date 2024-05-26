import editdistance
import numpy as np
import math
import time
from itertools import product
import csv
import matplotlib.pyplot as plt
import uuid
import os
import sys
os.chdir(sys.path[0])

from step_0_configure import num_trials, sd_min, sd_max, dwell_min, dwell_max, choice_list, dorado_call_auto, model, num_codebooks
from step_0_make_codebooks import codon_freq_dict

def load_list_of_lists_from_csv(filename):
    with open(filename, 'r') as csvfile:
        csv_reader = csv.reader(csvfile)
        data = [list(map(int, row)) for row in csv_reader]
    return data

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

# Retrieve the codebooks
codebooks = [None]*num_codebooks
for i in range(num_codebooks):
    with open(f'./1_codebooks/codebook_{i+1}.txt','r') as f:
        codebooks[i] = f.readlines()[0].split(' ')

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
        csv_filename = f'./3_squiggles/squiggles_sd{sd_mean_scaled}_K{dwell_mean}_cb{choice}.csv'
        all_signals = load_list_of_lists_from_csv(csv_filename)
        for i in range(num_trials):
            signal_int16 = np.array(all_signals[i])
            read_id = f"{str(sd_mean_scaled).zfill(8)}-{str(dwell_mean).zfill(4)}-{str(choice).zfill(4)}-0000-0000{str(i).zfill(8)}"
            start = time.time()
            rx_msg = scrappy.basecall_raw(signal_int16)[0]          # Normalisation of raw_data is irrelevant.'
            total_time += time.time()-start
            reads_dict[read_id][1] = rx_msg
    print(f'Scrappy basecalling finished, taking {total_time} seconds.')
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

stats = []

for params, edit_distances in results_dict.items():
    sd_mean = int(params[0:8]) / 100
    dwell_mean = int(params[9:13])
    choice = int(params[14:])
    
    metrics = evaluate_codebook(codebooks[choice-1])

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
    excel_summary = f"{choice},{len_msg_bases//3},{num_valid_trials},{sd_mean},{dwell_mean},{mean_result:.2f},{sd_result:.2f},{metrics[0]:.4f},{metrics[1]:.4f},{metrics[2]:.3f}"
    # print(f"Excel summary: {excel_summary}")
    summary_all.append(excel_summary)
    stats.append(list(metrics) + [Q_result])


# Output summary of trials
output_file = "summary.txt"
with open(output_file, "w") as f:
    for line in summary_all:
        f.write(line)
        f.write('\n')


stats = np.array(stats)
from sklearn.linear_model import LinearRegression

X = stats[:, :3]
y = stats[:, 3]

coefficients_list = []

for i in range(X.shape[1]):
    current_X = X[:, i].reshape(-1, 1)
    model = LinearRegression()
    model.fit(current_X, y)
    coefficients = (model.intercept_, model.coef_[0])
    coefficients_list.append(coefficients)

print(coefficients_list)

model = LinearRegression()
model.fit(X, y)

coefficients = model.coef_
intercept = model.intercept_

print("Coefficients:", coefficients)
print("Intercept:", intercept)

predicted_y = model.predict(X)

from sklearn.metrics import r2_score
# Assuming 'predicted_y' contains the predicted values and 'y' contains the actual values
r_squared = r2_score(y, predicted_y)
print("R-squared coefficient:", r_squared)

r_squared_list = []

# Calculate R-squared coefficient for each independent variable
for i in range(X.shape[1]):
    current_X = X[:, i].reshape(-1, 1)
    
    model = LinearRegression()
    model.fit(current_X, y)
    
    predicted_y = model.predict(current_X)
    
    r_squared = r2_score(y, predicted_y)
    r_squared_list.append(r_squared)

print("R-squared coefficients for each independent variable:", r_squared_list)

r_squared_list = []

# Calculate R-squared coefficient for each independent variable
for i in range(X.shape[1]):
    current_X = X[:, i].reshape(-1, 1)
    other_X = np.delete(X, i, axis=1)  # Remove the current independent variable
    
    other_X_cut = other_X[41:, :]
    other_y_cut = y[41:]
    model = LinearRegression()
    model.fit(other_X_cut, other_y_cut)
    
    predicted_y = model.predict(other_X)
    
    residuals = y - predicted_y
    
    r_squared = r2_score(y, y - residuals)  # R-squared using the residuals
    r_squared_list.append(r_squared)

print("R-squared coefficients for each independent variable with other impacts removed:", r_squared_list)