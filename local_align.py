from Matrix import *
from read_fasta import *

seq_list = read_fasta("protein-sequences.fasta")
seq1 = seq_list[0]
seq2 = seq_list[1]
slen1 = len(seq1)
slen2 = len(seq2)
mat1 = Matrix("BLOSUM50.txt")
print(seq1, seq2)

score_mat = []

for i in range(len(seq1) + 1):
    score_mat.append([0 for i in range(len(seq2) + 1)])


# create the direction mat
direc_mat = []
for i in range(len(seq1) + 1):
    direc_mat.append([[0, 0, 0] for i in range(len(seq2) + 1)])
# [0,1,0] means seq 2 has a gap
# [1,0,0] means seq1 has a gap
# [0,0,1] means diagonal
for i in range(1, len(seq2) + 1):
    direc_mat[0][i] = [0, 1, 0]
for i in range(1, len(seq1) + 1):
    direc_mat[i][0] = [1, 0, 0]
    # print(direc_mat)

def print_mat2(alist):
    for i in alist:
        for j in i:
            print('{0:4d}'.format(j), end=" ")
        print()


def print_mat3(alist):
    for i in alist:
        print(i)

def local_score(mat,gap):
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            t1 = score_mat[i - 1][j - 1] + mat.get(seq1[i - 1], seq2[j - 1])

            #[s(i,j-n(gap)+g(gap]
            t2 = []
            for k in range(j):
                t2.append(score_mat[i][k]+gap*(j-k))
            t3 =[]
            for k in range(i):
                t3.append(score_mat[k][j]+gap*(i-k))

            t2max = max(t2)
            t3max = max(t3)

            max_score= max(t1,t2max,t3max,0)
            score_mat[i][j] = max_score


local_score(mat1,-8)
print_mat2(score_mat)


