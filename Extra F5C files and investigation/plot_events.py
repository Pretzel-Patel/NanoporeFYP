
import os
import glob
import sys

os.chdir(sys.path[0])


import matplotlib.pyplot as plt
import pandas as pd

# Read the data from the text file
data = pd.read_csv('./5_aligned_events/events_reverse_1.tsv', sep='\t')


# Extract the first twenty rows

# Plotting
indices = [3,1,6,9, 10, 16, 17, 19, 22]

colours = ['m', 'r', 'b', 'g', 'skyblue']

N = 2
fig, axs = plt.subplots(N,N)
for i in range(N**2):
    data_sub = data[data.read_index == indices[i]]


    # Plotting the rectangular signal
    curr_time = 0
    prev_pos = 0
    ax = axs[i//N,i%N]
    for index, row in data_sub.iterrows():
        event_level_mean = row['event_level_mean']
        event_length = row['event_length']
        pos = row['position']
        if pos != prev_pos:
            prev_pos = pos
            ax.plot([curr_time, curr_time], [40, 140], 'k--')
        ax.plot([curr_time, curr_time+event_length], [event_level_mean, event_level_mean], color=colours[pos % len(colours)])
        curr_time += event_length
        if pos >= 10:
            break

    ax.plot([curr_time, curr_time], [40, 140], 'k--')
    # Set labels and title
    ax.set_xlabel('Time')
    ax.set_ylabel('Current Level')
    ax.set_title(f'Read {indices[i]}')

plt.show()