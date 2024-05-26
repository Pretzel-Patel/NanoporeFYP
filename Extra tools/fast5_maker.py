'''
fast5_maker.py

This file takes a list of current signals, and stores this information into a FAST5 file.

It also shows all the metadata needed to create a FAST5 file.
This requires installing the fast5_research package, which worked on python 3.6 in Linux.
'''

from fast5_research import Fast5
import matplotlib.pyplot as plt
import numpy as np
import uuid
import os
import sys
os.chdir(sys.path[0])

# Conversion coefficients used for digitisation
# current_pA = scale * (raw_int + offset) = sd_current * current_normalised + mean_current = sigma * current_normalised + mu
mu = 66.27
sigma = 12.19

offset = 10
scale = 0.1755


# Create 10 random squiggles with random read_ids
current_pA = np.random.uniform(50, 140, 4000)

# Example of how to digitize current data
start, stop = int(min(current_pA - 1)), int(max(current_pA + 1))
rng = stop - start
digitisation = 8192.0
bins = np.arange(start, stop, rng / digitisation)
# np.int16 is required, the library will refuse to write anything else
signal = np.digitize(current_pA, bins).astype(np.int16)

# Name of new FAST5 file
fast5_filename= 'new.fast5'

# The following are required meta data
channel_id = {
    'digitisation': digitisation,
    'offset': 0,
    'range': rng,
    'sampling_rate': 4000,
    'channel_number': 1,
    }
read_id = {
    'start_time': 0,
    'duration': len(signal),
    'read_number': 1,
    'start_mux': 1,
    'read_id': str(uuid.uuid4()),           # Generate a random UUID
    'scaling_used': 1,
    'median_before': 0,
}
tracking_id = {
    'exp_start_time': '1970-01-01T00:00:00Z',
    'run_id': str(uuid.uuid4()).replace('-',''),
    'flow_cell_id': 'FAH00000',
}
context_tags = {}

with Fast5.New(fast5_filename, 'w', tracking_id=tracking_id, context_tags=context_tags, channel_id=channel_id) as h:
    # print(raw_data)
    h.set_raw(signal, meta=read_id, read_number=1)
