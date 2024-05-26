import editdistance


w1 = 'ATGTGTCCTGTATTTGGTTCAGTTGCGTCTTGTCGGTGCTAGTATTGGTTGGCTTGGCTGTGGGCATCATTTCCATAGCCTAGCATCATTTCCCGTTATAGCGGCCTAGTCCACATTTCCCAGCATCATTTGCCAGCCGTTGAGCATCATTTCCAAGACTCCCTCCAGTCATTAAAGACTCCCACCTGTTTTTTGGTTA'
w2 = 'TGTTCGGTGCTCGTATTCGTAGGCTTGGCTGTGGGCATCATTTCCCCATAGCCTAGCATCATTTCCCGTTATAGCAGCATCATTTCCACATTTCCCAGCATCATTTCCCAGCCGTTGAGCATCATTTCCAAGACTCCCTCCAGTCATTAAAGACTCCGCCTGTTTTGGTTAAACACCCAAGC'
# w1 = 'TACGTATGCTTCGTATTACT'
# w2 = 'GCTTGGGTGTTTAACCAAAA'
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



print(f"My code found edit distance of {a[m][n][0]}. [S, D, I] = {a[m][n][1::]}")
print(f'Package found edit distance of {editdistance.eval(w1, w2)}.')