# Empty the folders
import os
import sys
import glob

os.chdir(sys.path[0])

folders = ['2_tx_msg', '3_squiggles', '4_reads', '5_rx_msg']
for folder in folders:
    files = glob.glob(f'{folder}/*')
    for f in files:
        os.remove(f)