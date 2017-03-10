# in a non recursive way
from Matrix import *
from read_fasta import *

seq_list = read_fasta("WW-sequence.fasta")
seq1 = seq_list[0][:10]
seq2 = seq_list[1][:10]
slen1 = len(seq1)
slen2 = len(seq2)
mat1 = Matrix("BLOSUM50.txt")
print(seq1, seq2)

score_mat = []

for i in range(len(seq1) + 1):
    score_mat.append(['#' for i in range(len(seq2) + 1)])
score_mat[0][0] = 0
t = 0
for i in range(1, len(seq2) + 1):
    t += mat1.get('*', seq2[i - 1])
    score_mat[0][i] = t
t = 0
for i in range(1, len(seq1) + 1):
    t += mat1.get('*', seq1[i - 1])
    score_mat[i][0] = t

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


def score2(mat):

    for i in range(1,len(seq1)+1):
        for j in range(1,len(seq2)+1):

            t1 = score_mat[i-1][j]-5

            t2 = score_mat[i][j-1]-5

            t3 = score_mat[i-1][j-1] + mat.get(seq1[i-1],seq2[j-1])

            maxscore = max(t1, t2, t3)
            score_mat[i][j] = maxscore
            if (t1 == maxscore):
                direc_mat[i][j][1] = 1
            if (t2 == maxscore):
                direc_mat[i][j][0] = 1
            if (t3 == maxscore):
                direc_mat[i][j][2] = 1
    return score_mat[len(seq1)][len(seq2)]


print(score2(mat1))
def print_mat2(alist):
    for i in alist:
        for j in i:
            print('{0:4d}'.format(j),end=" ")
        print()

def print_mat3(alist):
    for i in alist:
        print(i)
print_mat2(score_mat)
print_mat3(direc_mat)
