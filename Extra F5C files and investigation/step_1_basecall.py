'''



'''
import os
import sys


# dorado basecaller -v --emit-fastq --batchsize 64 dna_r9.4.1_e8_fast@v3.4 4_reads/all_reads.pod5 > 5_rx_msg/all_calls.fa
from step_0_configure import num_files, model, dorado_call_auto, direction

os.chdir(sys.path[0])

for n in range(1,9):
    sequence = f'{direction}_{n}'
    # Check if basecalled file exists. If not, then basecall the fasta file to fastq.
    command = f'dorado basecaller -v --emit-fastq --batchsize 64 --read-ids ./2_read_ids/{sequence}.txt ../models_dorado/{model} ./1_input_fast5 > ./3_basecalls_fa/{sequence}.fa'
    # Run dorado in command line using python
    if dorado_call_auto:
        os.system(command)
    else:
        print(f'Run this in windows shell\n{command}')

    # Check if basecalled file exists. If not, then basecall the fasta file to bam.
    command = f'dorado basecaller -v --batchsize 64 --read-ids ./2_read_ids/{sequence}.txt --reference ./2_read_ids/reference_{sequence}.fa ../models_dorado/{model} ./1_input_fast5 > ./4_basecalls_bam/{sequence}.bam'
    # Run dorado in command line using python
    if dorado_call_auto:
        os.system(command)
    else:
        print(f'Run this in windows shell\n{command}')
