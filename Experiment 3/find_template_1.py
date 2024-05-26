import pandas as pd
from ont_fast5_api.fast5_interface import get_fast5_file, fast5_subset
import time
import os
import sys
import editdistance


os.chdir(sys.path[0])


df = pd.read_csv('4_basecalls_csv/all_reads_full_simple.csv')
df.set_index('read_id', inplace=True)

fast5_filepath = './1_input_fast5/240129_R10_newsystem_verified0.fast5'

def print_all_raw_data():
    with get_fast5_file(fast5_filepath, mode="r") as f5:
        for read in f5.get_reads():
            if (df.loc[read.read_id, 'min_col']) == 'r1dist':
                raw_data = read.get_raw_data()
                print(read.read_id, raw_data)


print_all_raw_data()