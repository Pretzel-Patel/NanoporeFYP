"""
pod5_peeker.py 

This file peeks into a POD5 file and allows the user to inspect its data.

By default, the first 10 read_ids will be printed. This can be switched off using the variable parameter.
It is recommended that a debugger is used to examine the content of the 'read' variable.
This requires installing the POD5 package, which works on windows.
"""

import matplotlib.pyplot as plt
import numpy as np
import pod5 as p5

import os
import sys

os.chdir(sys.path[0])

# Using the sample.pod5 file provided
pod5_file = "sample.pod5"
selected_read_id = "00000010-0010-0002-0000-000000000009"
print_read_ids = True

with p5.Reader(pod5_file) as reader:
    # Print first 10 read_ids
    if print_read_ids:
        print('Here are the first 10 read_ids found in the file.')
        for read_id in reader.read_ids[:10]:
            print(read_id)
        print()
    
    # Read the selected read from the pod5 file
    # next() is required here as Reader.reads() returns a Generator
    read = next(reader.reads([selected_read_id]))

    # Get the signal data and some other run info.
    sample_rate = read.run_info.sample_rate
    offset = read.calibration.offset
    scale = read.calibration.scale
    signal = read.signal
    signal_pa = read.signal_pa

    ## INSERT DEBUG POINT HERE AND EXAMINE 'read' VARIABLE

    # Print signal statistics for digitised and original signal
    print(f'Digitised signal: Min: {min(signal)}, Mean: {signal.mean():.2f}, Max: {max(signal)}.')
    print(f'Current signal (pA): Min: {min(signal_pa):.2f}, Mean: {signal_pa.mean():.2f}, Max: {max(signal_pa):.2f}.')

    # Plot signal
    plt.plot(signal_pa)
    plt.xlabel("Sample number")
    plt.ylabel('Current level (pA)')
    plt.title('Current value vs. sample number')
    # plt.savefig("pod5data.png")
    plt.show()
    
