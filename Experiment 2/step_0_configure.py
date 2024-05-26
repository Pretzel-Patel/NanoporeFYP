# Trial constants
len_msg = 200                       # number of codewords
PAD_str = 'ACGTA'
num_trials = 50                   # number of reads to generate for each combination of parameters

# Trial parameters
sd_min = 10                         # 10 means 0.10 normalised noise
sd_max = 10
dwell_min = 10                      # dwell time in number of bases
dwell_max = 10

pop_codeboooks = 2000
GC_codebooks  = 0
homo_codebooks = 0
num_codebooks = pop_codeboooks+GC_codebooks+homo_codebooks
choice_list = list(range(1,num_codebooks+1))             # index of codebooks to test, corresponding to numbering in ./1_codebooks/


# Squiggle generation method
squiggle_generation = 'kmer'        # 'kmer' | 'scrappie'

# Allow files to automatically call each other and run in command line
dorado_call_auto = True

# Basecalling model to be used.
# R9.4 Scrappie
# model, kmer_filename, sample_rate =  'scrappie',      None,   4000

# R9.4, Dorado
# model, kmer_filename, sample_rate = 'dna_r9.4.1_e8_fast@v3.4',     'R94_450bps_6mer.model',    4000
# model, kmer_filename, sample_rate = 'dna_r9.4.1_e8_hac@v3.3',      'R94_450bps_6mer.model',    4000
# model, kmer_filename, sample_rate = 'dna_r9.4.1_e8_sup@v3.3',      'R94_450bps_6mer.model',    4000

# R10.4, 400bps, 4.2.0
# model, kmer_filename, sample_rate = 'dna_r10.4.1_e8.2_400bps_fast@v4.2.0',     'R104_400bps_9mer.model',   5000
# model, kmer_filename, sample_rate = 'dna_r10.4.1_e8.2_400bps_hac@v4.2.0',      'R104_400bps_9mer.model',   5000
model, kmer_filename, sample_rate = 'dna_r10.4.1_e8.2_400bps_sup@v4.2.0',      'R104_400bps_9mer.model',   5000

# R10.4, 260bps, 4.1.0
# model, kmer_filename, sample_rate = 'dna_r10.4.1_e8.2_260bps_fast@v4.1.0',     'R104_260bps_9mer.model',   4000
# model, kmer_filename, sample_rate = 'dna_r10.4.1_e8.2_260bps_hac@v4.1.0',      'R104_260bps_9mer.model',   4000
# model, kmer_filename, sample_rate = 'dna_r10.4.1_e8.2_260bps_sup@v4.1.0',      'R104_260bps_9mer.model',   4000

# kmer squiggle generation setup
if squiggle_generation == 'scrappie':
    # No kmer parameters needed for scrappie
    kmer_len = None
    kmer_centre = None
    kmer_legacy = None
elif kmer_filename.split('_')[0] == 'R104':
    kmer_len = 9
    kmer_centre = 5                 # For 9mer, this is 5 because the 9mer level mean corresponds to the 5th base in the 9mer.
    kmer_legacy = False             # kmer legacy files have 6 columns, while new kmer files only have 2 columns
elif kmer_filename.split('_')[0] == 'R94':
    kmer_len = 6
    kmer_centre = 3                 # For 6mer, this is 3 because the 6mer level mean corresponds to the 3rd base in the 6mer.
    kmer_legacy = True              # kmer legacy files have 6 columns, while new kmer files only have 2 columns
else:
    raise Exception('Error in param setup')



# Offset pA values (for R94) - Not important
mu = 66.27
sigma = 12.19

# Offset for int16 values -  Not important
offset = 10
scale = 0.1755