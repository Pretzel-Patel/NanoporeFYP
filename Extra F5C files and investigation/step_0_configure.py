# Input FASTA specification
input_fasta_prefix = "240129_R10_newsystem_verified"
num_files = 52           # 52
direction = 'reverse'

read_length  = 182

# Offset pA values (for R94)
mu = 66.27
sigma = 12.19

# Offset for int16 values
offset = 10
scale = 0.1755

# Adrians function
def reverse_complement(s):
    d = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    return ''.join([d[c] for c in reversed(s)])

# Trying to improve basecalling
ignore_bases = 2
ont_seq = 'GCTTGGGTGTTTAACCAAAA'
ont_seq_rev = 'TTTTGGTTAAACACCCAAGC' # reverse_complement(ont_seq)

# Model to be used. Need to change kmer setup as well. and possibly kmer_legacy, kmer_len, kmer_centre
# R9.4
# model, kmer_filename, sample_rate = 'dna_r9.4.1_e8_fast@v3.4',     'R94_450bps_6mer.model',    4000
# model, kmer_filename, sample_rate = 'dna_r9.4.1_e8_hac@v3.3',      'R94_450bps_6mer.model',    4000
# model, kmer_filename, sample_rate = 'dna_r9.4.1_e8_sup@v3.3',      'R94_450bps_6mer.model',    4000
# R10.4, 400bps, 4.2.0
# model, kmer_filename, sample_rate = 'dna_r10.4.1_e8.2_400bps_fast@v4.2.0',     'R104_400bps_9mer.model',   5000
# model, kmer_filename, sample_rate = 'dna_r10.4.1_e8.2_400bps_hac@v4.2.0',      'R104_400bps_9mer.model',   5000
model, kmer_filename, sample_rate = 'dna_r10.4.1_e8.2_400bps_sup@v4.2.0',      'R104_400bps_9mer.model',   5000
# R10.4, 400bps, 4.2.0
# model, kmer_filename, sample_rate = 'dna_r10.4.1_e8.2_260bps_fast@v4.1.0',     'R104_260bps_9mer.model',   4000
# model, kmer_filename, sample_rate = 'dna_r10.4.1_e8.2_260bps_hac@v4.1.0',      'R104_260bps_9mer.model',   4000
# model, kmer_filename, sample_rate = 'dna_r10.4.1_e8.2_260bps_sup@v4.1.0',      'R104_260bps_9mer.model',   4000

# Using manual dorado call in command line?
dorado_call_auto = True
# kmer setup
if kmer_filename.split('_')[0] == 'R104':
    kmer_len = 9
    kmer_centre = 5         # For 9mer, this is 5 because the 9mer level mean corresponds to the 5th base in the 9mer.
    kmer_legacy = False   # kmer legacy have 6 columns, while other kmer only have 2 columns
elif kmer_filename.split('_')[0] == 'R94':
    kmer_len = 6
    kmer_centre = 3         # For 6mer, this is 3 because the 6mer level mean corresponds to the 3rd base in the 6mer.
    kmer_legacy = True   # kmer legacy have 6 columns, while other kmer only have 2 columns
else:
    raise Exception('Error in param setup')
