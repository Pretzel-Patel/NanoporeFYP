'''
This file should be run on a linux machine which has samtools and f5c installed.
The GPU machine in room G16 encounters a stack smash error when f5c is used. 
This uses the default minimum mapping quality of 20.
There is an option to print read_ids instead of index
See here: https://hasindu2008.github.io/f5c/docs/output#eventalign
    for explanation of output
'''
import os
import sys

from step_0_configure import direction

# dorado basecaller -v --emit-fastq --batchsize 64 dna_r9.4.1_e8_fast@v3.4 4_reads/all_reads.pod5 > 5_rx_msg/all_calls.fa
# from step_0_configure import 
os.chdir(sys.path[0])

for n in range(1,9):
    sequence = f'{direction}_{n}'
    os.system(f'samtools sort ./4_basecalls_bam/{sequence}.bam > ./4_basecalls_bam/{sequence}_sorted.bam')
    os.system(f'samtools index ./4_basecalls_bam/{sequence}_sorted.bam')
    os.system(f'~/f5c-v1.4/f5c index -d 1_input_fast5/ ./3_basecalls_fa/{sequence}.fa')
    os.system(f'~/f5c-v1.4/f5c eventalign --samples -b ./4_basecalls_bam/{sequence}_sorted.bam -g reference.fa -r ./3_basecalls_fa/{sequence}.fa --pore r10 > ./5_aligned_events/events_{sequence}.tsv')
    
