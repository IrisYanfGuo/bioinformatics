# in a non recursive way
from Matrix import *
from read_fasta import *

# affine gap penalty
affine = -0
gap = -8

seq_list = read_fasta("protein-sequences.fasta")
seq1 ='THISLINE'
seq2 = 'ISALIGNED'
slen1 = len(seq1)
slen2 = len(seq2)
mat1 = Matrix("BLOSUM62.txt")
print(seq1)
print(seq2)

score_mat = []
direc_mat = []
for i in range(len(seq1) + 1):
    direc_mat.append([[0, 0, 0,0] for i in range(len(seq2) + 1)])

for i in range(len(seq1) + 1):
    score_mat.append([0 for i in range(len(seq2) + 1)])
# initialize


def print_mat2(alist):
    for i in alist:
        for j in i:
            print('{0:4d}'.format(j), end=" ")
        print()


def print_mat3(alist):
    for i in alist:
        print(i)


def local_score(mat):
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            t1 = score_mat[i - 1][j - 1] + mat.get(seq1[i - 1], seq2[j - 1])

            # [s(i,j-n(gap)+g(gap]
            t2 = []
            for k in range(j - 1, -1, -1):
                t = score_mat[i][k] + gap * (j - k) + int((j-k)*(j-k-1)*affine/2)
                if (t < 0):
                    break
                t2.append(t)

            t3 = []
            for k in range(i - 1, -1, -1):

                t = score_mat[k][j] + gap * (i - k) + int((j-k-1)*(j-k)*affine/2)
                if (t < 0):
                    break
                t3.append(t)

            if len(t2) > 0:
                t2max = max(t2)
            else:
                t2max = -1
            if len(t3) > 0:
                t3max = max(t3)
            else:
                t3max = -1

            max_score = max(t1, t2max, t3max, 0)
            score_mat[i][j] = max_score

            if max_score > 0:
                if max_score == t1:
                    direc_mat[i][j][2] = 1

                elif max_score == t2max:
                    # horizental
                    direc_mat[i][j][1]=1
                else:
                    # vertical, seq1 has a gap
                    direc_mat[i][j][0]=1

a,i = max([1,2,3]),argmax([1,2,3])




local_score(mat1)
print_mat2(score_mat)
print_mat3(direc_mat)

