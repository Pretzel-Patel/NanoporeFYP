'''
This iteration deletes bases from the dorado reads when the ont_seq or ont_seq_rev is detected
ont_seq belongs to template sequences, and I note that around 30 bases are deleted on average.
ont_seq_rev belongs to reverse sequences, and I note that around 10 bases are deleted on average.
Some further investigation could be done plotting s1 vs. score1 and e2 vs score1.
    This could give an indication of how many bases should be deleted.
    
'''
import os
import sys
import pandas as pd
import editdistance


os.chdir(sys.path[0])
# dorado basecaller -v --emit-fastq --batchsize 64 dna_r9.4.1_e8_fast@v3.4 4_reads/all_reads.pod5 > 5_rx_msg/all_calls.fa
from step_0_configure import num_files, model, dorado_call_auto, input_fasta_prefix, ont_seq, ont_seq_rev
from my_lookup import search


def my_search(word1, word2, where):
    '''
    Look for the occurence of word1 in word2.
    Return the location of the occurence and the editdistance.
    my_search is much faster than search
    ''' 
    # Asymetric as start of read is more uncertain
    if where == 'start':
        return search(word1, word2[:60])
    elif where == 'end':
        N = len(word2)
        (s, e, score) = search(word1, word2[-45:])
        return (s+N-45, e+N-45, score)
    else:
        raise ValueError


# Check if basecalled file exists. If not, then basecall the fasta file.
for file_num in range(num_files):
    fast5_filename = f"./1_input_fast5/{input_fasta_prefix}{file_num}.fast5"
    if not os.path.isfile(fast5_filename):
        command = f'dorado basecaller -v --emit-fastq --batchsize 64 ../models_dorado/{model} {fast5_filename} > ./3_basecalls_fa/calls_{file_num}.fa'
        # Run dorado in command line using python
        if dorado_call_auto:
            os.system(command)
        else:
            print(f'Run this in windows shell\n{command}')




# Assuming your table is stored as a dictionary where keys are read_id and values are DNA sequences
# File names
reads_dict = dict()
k = 0
for file_num in range(num_files):
    rx_msg_file = f'./3_basecalls_fa/calls_{file_num}.fa'

    if dorado_call_auto:
        with open(rx_msg_file,'r') as f:
            rx_msg_all = f.readlines()
    else:
        with open(rx_msg_file,'r',encoding='utf-16le') as f:       # Use this when command line dorado is used
            rx_msg_all = f.readlines()

    i = 0
    # rx_msg_all = [line.replace('\x00','') for line in rx_msg_all]
    while i < len(rx_msg_all)-1:
        # print(i)
        assert '@' in rx_msg_all[i]
        read_id = (rx_msg_all[i].split('@'))[1][:-1]
        
        if reads_dict.get(read_id) != None:
            print(f'Double ID: {read_id} in file {file_num}')
            sys.exit()
        sequence = rx_msg_all[i+1][:-1]
        # Remove ont_seq or its reverse from sequence
        (s1, e1, score1) = my_search(ont_seq, sequence, 'start')
        (s2, e2, score2) = my_search(ont_seq_rev, sequence, 'end')
        score = min(score1, score2)
        # A problem is that a match cannot be found in every strand, and we cannot even know which side to delete bases from.
        # Therefore, I propose that we do not delete the ont_seq itself - only the stuff preceding or succeeding it.
        # print(i//4)
        if score <= 3 and score1 != score2:
            if score1 == score:
                N = s1
                sequence = sequence[s1:]
                print(f'{k} template detected: {s1}:{e1}. {N} bases deleted with score {score1}')
            else:
                N = len(sequence) - e2
                sequence = sequence[:e2]
                print(f'{k} reverse  detected: {s2}:{e2}. {N} bases deleted with score {score2}')
        reads_dict[read_id] = sequence
        i += 4
        k += 1

# Import reference sequences
with open('./ref_seq.txt','r') as f:
    rows = f.readlines()
ref_list = []
for i in range(1, len(rows)):
    row = rows[i].strip().split('\t')
    # Remove ont_seq or its reverse. Presumably no reference sequence has both of these. Changed my mind and am no longer doing this.
    # row[2] = row[2].replace(ont_seq, '')
    # row[2] = row[2].replace(ont_seq_rev, '')
    ref_list.append(row)
ref_seq = pd.DataFrame(ref_list, columns=['orientation', 'label', 'sequence'])



# Store as CSV file
print(f'Total number of reads is {len(reads_dict)}')
df = pd.DataFrame(list(reads_dict.items()), columns=['read_id', 'DNA_sequence'])
for i in range(len(ref_seq)):
    compare = f'{ref_seq.orientation.loc[i][0]}{ref_seq.label.loc[i]}dist'
    print(compare)
    df[compare] = df.apply(lambda x: editdistance.eval(x['DNA_sequence'],ref_seq.sequence.loc[i]), axis=1)
df['min_col']=df.loc[:,'t1dist':'r8dist'].idxmin(1)
df['orientation'] = df.apply(lambda x: 'template' if x['min_col'][0] == 't' else 'reverse', axis=1)
df['label'] = df.apply(lambda x: int(x['min_col'][1]), axis=1)
df['min_1'] = df.loc[:,'t1dist':'r8dist'].min(axis=1)
df['min_2'] = second_smallest_values = df.loc[:,'t1dist':'r8dist'].apply(lambda row: sorted(row)[1], axis=1)
df['strength'] = (df['min_2']/df['min_1']).round(4)

output_df = df[['read_id', 'orientation', 'label', 'strength', 'DNA_sequence']]
output_df.to_csv(f'4_basecalls_csv/all_reads_summary.csv')
df.to_csv(f'4_basecalls_csv/all_reads_full_ont_removed.csv')