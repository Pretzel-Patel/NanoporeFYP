This pipeline was largely developed in late May 2024.
It is a combination of previous pipelines, which correspond to experiment 1

Its purpose is to take a selection of codebooks, a range of sd, a range of dwell, and N, a number of trials, and then go through each combination and create N different tx_messages, squiggles and basecalls using dorado. 
The basecalls are single read, and are compared to the tx_message to determine the error rate and Q-score. 

Squiggle generation is via Scrappie or k-mer tables.

Basecalling is via Scrappie or Dorado

Preparation for running trials:
1. Ensure models_dorado and models_kmer is in the same directory as the Experiment 1 folder.
2. Adjust parameters in the step_0_configure file are correct.

Steps to running trials:
1.  Run step_0_clear.py in python.
2.  Run step_1_make_squiggle.py. If dorado_call_auto is True, then step_2 and step_3 will run automatically.
3.  Run step_2_make_pod_5_reads.py
3A. If dorado_call_auto is True, then this will call dorado in-file, and also call step_3, so you are done.
3B. If dorado_call_auto is False, then this will print a dorado call, which should be run in windows.
4.  Run step_3_summarise

Notes:
If Scrappie is used for squiggle generation, then dorado_call_auto should be False and step_1 should be run in Linux.
If Scrappie is used for basecalling, then step_2 and step_3 should be run in Linux.