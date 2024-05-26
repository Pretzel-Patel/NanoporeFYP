'''
This file contains a function which looks for the best substring match of a word1 in a (presumably longer) word2.
The location and editdistance of the match are returned.
The purpose is to use this to find and remove the ont_seq from the dorado reads, so that it is easier to find the reference sequence.
It is expected that template and reverse reads are less likely to be mixed up.
Comparing two full-length (200base) sequences takes 2 seconds. 
'''


import editdistance

import time
import os
import sys
import pandas as pd
import editdistance


os.chdir(sys.path[0])


def search(word1, word2):
    '''
    Look for the occurence of word1 in word2.
    Return the location of the occurence and the editdistance.
    ''' 
    best_start = 0
    best_end = 0
    best_count_overall = len(word2)+len(word1)
    for start_candidate in range(len(word2)):
        w1 = word1
        w2 = word2[start_candidate:]
        len1 = len(w1)
        len2 = len(w2)
        a =[[[0,0,0,0]]*(len2+1) for _ in range(len1+1)]

        for m in range(0,len1+1):
            a[m][0]=[m,0,m,0]       # deletion errors
            
        for n in range(0,len2+1):
            a[0][n]=[n,0,0,n]       # insertion errors

        ins_count = 0
        del_count = 0
        sub_count = 0  

        best_count = len1+len2
        good_end = 0
        for m in range (1,len1+1):
            for n in range(1,len2+1):
                if w1[m-1]==w2[n-1]:
                    a[m][n] = a[m-1][n-1]
                else :
                    ins_score = a[m][n-1][0]
                    del_score = a[m-1][n][0]
                    sub_score = a[m-1][n-1][0]
                    if ins_score < del_score and ins_score < sub_score:
                        # ins_count += 1
                        a[m][n] = a[m][n-1][:]
                        a[m][n][0] += 1
                        a[m][n][3] += 1
                    elif del_score < sub_score:
                        # del_count += 1
                        a[m][n] = a[m-1][n][:]
                        a[m][n][0] += 1
                        a[m][n][2] += 1
                    else:
                        # sub_count += 1
                        a[m][n] = a[m-1][n-1][:]
                        a[m][n][0] += 1
                        a[m][n][1] += 1
                if m == len1 and a[m][n][0] < best_count:
                    best_count = a[m][n][0]
                    good_end = n
                    # I don't care where it starts, I care where the match ends.
        if best_count < best_count_overall:
            best_count_overall = best_count
            best_start = start_candidate
            best_end = best_start + good_end
    return (best_start, best_end, best_count_overall)




if __name__ == '__main__':
    # Import reference sequences
    with open('./ref_seq.txt','r') as f:
        rows = f.readlines()
    ref_dict = {}
    for i in range(1, len(rows)):
        row = rows[i].strip().split('\t')
        # Remove ont_seq or its reverse. Presumably no reference sequence has both of these. Changed my mind and am no longer doing this.
        # row[2] = row[2].replace(ont_seq, '')
        # row[2] = row[2].replace(ont_seq_rev, '')
        ref_dict[row[0]+row[1]] = row[2]
    # ref_seq = pd.DataFrame(ref_list, columns=['orientation', 'label', 'sequence'])
        
    start_time = time.time()
    results_dict = {}
    for key, value in ref_dict.items():

        word2 = 'ATGTACTCGTTAGTTTACGTATTGCTTGTTTAGGTGCTAGTATTCGTAGGCTTGGCTGTGGGCATCGTTCCCATAGCCTAGCATCATTTGCCCGTTATAGCGGCGTCATTCCCACATTTTCCAGCATCATTTCCCAGCCGTTGAGCATAGTTTCCAAGACTCCCTCCAGTCATTAAAGACTCCGCCTGTTTTGGTTCAACGCCCGAGCAGCAATACGTGG'
        word2 = 'GTTCCGTATTACTCGGTGGGTCGATAGGCTTCAGCTATGGGCATCATTTCCACAGATAGCAGCATCATTTCCCATAAAGCCAGCATCATTTCCGTTAGTGAGAGCATCATTTCCAATTCCGAGAGCATGATTTCCAAGACTCCCTCCAGTCACAGAGCCGTTGAGAGATTTTCGTTAT'
        word2 = 'CCCATTCCGATAGCATCATTTCCTTTCAGAGCAGCATCATTTCCCATAGCCACAGCATCATTTCCGTTTCCGATAGCATCCAAGACTCCCTCCAGTCAGAAGGATCCAT'
        word1 = value

        results_dict[key] = search(word1, word2)
    # print(f'{time.time()-start_time}')
    # # print(f"My code found edit distance of {a[m][n][0]}. [S, D, I] = {a[m][n][1::]}")
    # # print(f'Package found edit distance of {editdistance.eval(w1, w2)}.')
    # print(f'Cutoffs identified as {best_start}:-{len(word2)-best_end}, and best score as {best_count_overall}.')
    # print(f"Matched '{word1}' with '{word2[best_start:best_end]}' from '{word2}'.")
    print(f'Start (0), end ({len(word2)-1}), editdistance')
    for key, value in results_dict.items():
        print(f'{key}: {value}')
