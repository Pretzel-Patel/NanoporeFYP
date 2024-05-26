"""
fast5_peeker.py 

This file peeks into a FAST5 file and allows the user to inspect its data.

Full documentation is available on https://nanoporetech.github.io/fast5_research/index.html
This requires installing the fast5_research package, which worked on python 3.6 in Linux.
"""

from fast5_research import Fast5
import matplotlib.pyplot as plt
import numpy as np
import uuid
import os
import sys
os.chdir(sys.path[0])

# Enter filename
filename = 'sample.fast5'

with Fast5(filename) as fh:
    raw = fh.get_read(raw=True)
    # fa = fh.get_reference_fasta()
    summary = fh.summary()

print('Raw is {} samples long.'.format(len(raw)))
print(f'Mean is {np.mean(raw)} and std is {np.std(raw)}')
print(raw)
# print('Summary {}.'.format(summary))
