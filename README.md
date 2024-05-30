# FYP on Coding Problems in Nanopore Sequencing

This repository contains the code used to perform each of the experiments in my Final Year Project at Monash University.

The _Report and Presentation_ directory contains my Final Report which investigates the use of codons to encode data before transmission through a DNA storage channel with nanopore sequencing. This report contains:
* A description of the method and data collection process.
* An overview of FASTA, FAST5 and POD5 file formats as an appendix.
* Guides to using Scrappie, Docker and F5C as appendices.

This repository contains the following code:
* Experiment 1 - A simulation pipeline which can be used to test codebooks against a range of padding, noise, dwell times, squiggle-generation methods, basecallers and basecalling models.
* Experiment 2 - A simulation pipeline which can be used to create and test a large number of codebooks, and observe the relationship between the quality score and three codebook metrics.
* Experiment 3 - A classifier which takes in a large number of reads and classifies these as one of 16 reference sequences.
* Extra F5C files and investigation - Contains several python notebooks explaining the F5C event alignment process and how this can be used to generate a consensus squiggle.
* Extra tools - Python files to read and generate FAST5 and POD5 files.
* models_dorado - A directory to hold dorado basecalling models.
* models_kmer - A directory to hold k-mer tables.
