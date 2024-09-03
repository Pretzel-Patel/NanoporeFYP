# Make kmer dictionary which generalises what scrappie does.
# .model files should not have headers.

from step_0_configure import kmer_filename, kmer_legacy

import os
import sys

os.chdir(sys.path[0])

kmer_dict = {}

if kmer_legacy:
    mean_pA = 90.5 # In kmer table, level_mean = 90.5 corresponds to scrappie mean 0, 
    scale_pA = 14.44 # in kmer table, level_mean = 61.5 corresponds to scrappie mean -2, and level_mean = 119.7 corresponds to scrappie mean +2.

else:
    mean_pA = 0 # No scaling needed for new tables
    scale_pA = 1

with open(f"../models_kmer/{kmer_filename}") as fileobject:
    for line in fileobject:
        line = line.split('\t')
        level_mean = float(line[1])
        # level_std = float(line[2])
        # sd_mean = float(line[3]) # not used
        # sd_std = float(line[4]) # not used
        kmer_dict[line[0]] = (level_mean-mean_pA)/scale_pA

