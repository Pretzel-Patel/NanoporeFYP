"""
Example use of pod5 to plot the signal data from a selected read.
can also pause debugger and look inside read variable
"""

import matplotlib.pyplot as plt
import numpy as np
import pod5 as p5
from step_0_configure import file_num, model


# Import pod5 file
filename = f"AV_{file_num}"
pod5_filename = f"2_converted_pod5/{filename}.pod5"
selected_read_id_0 = "00000016-0012-0001-0000-000000000183"
selected_read_id = selected_read_id_0


# 10 Nov tests multi
# trial = int(input("Enter trial: "))
# example_pod5 = f"10nov_testing_pod5/multi_{trial}.pod5"
# selected_read_id = "00000000-0000-0000-0000-000000000001"

with p5.Reader(pod5_filename) as reader:

    
    # Read the selected read from the pod5 file
    # next() is required here as Reader.reads() returns a Generator
    read = next(reader.reads([selected_read_id]))

    # Get the signal data and sample rate
    # sample_rate = read.run_info.sample_rate
    # signal = read.signal
    signal = list(read.signal_pa)
    print(type(read.run_info.protocol_start_time))
    # # Compute the time steps over the sampling period
    # time = np.arange(len(signal)) / sample_rate
    # import csv
    # with open('scrappie/testing/single_data.csv', 'w', newline='') as file:
    #     writer = csv.writer(file)
    #     writer.writerow(list(read.signal))
    # plt.plot(time, signal)
    # plt.show()
    print(np.mean(signal))
    print(min(signal))
    print(max(signal))
    plt.plot(list(range(len(signal))),signal)
    plt.xlabel("Samples")
    plt.ylabel('int current value')
    plt.title('Current value vs. sample number')
    plt.savefig("pod5data.png")
    plt.show()
    offset = read.calibration.offset
    scale = read.calibration.scale
    signal_pA = []
    for i in range(len(signal)):
        signal_pA.append((signal[i]+offset)*scale)
    
    pass


