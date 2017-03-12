# in a non recursive way
from Matrix import *
from read_fasta import *

# affine gap penalty
E = -4
I = -6 # affine gap


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
    direc_mat.append([[0, 0, 0,I,I] for i in range(len(seq2) + 1)])

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


def local_score(mat,istart=1,jstart=1):
    for i in range(istart, len(seq1) + 1):
        for j in range(jstart, len(seq2) + 1):
            t1 = score_mat[i - 1][j - 1] + mat.get(seq1[i - 1], seq2[j - 1])

            # [s(i,j-n(gap)+g(gap]

            t2 =score_mat[i][j-1] +direc_mat[i][j-1][4]



            t3 = score_mat[i-1][j] + direc_mat[i-1][j][3]



            max_score = max(t1, t2, t3, 0)
            score_mat[i][j] = max_score

            if max_score > 0:
                if max_score == t1:
                    direc_mat[i][j][2] = 1
                    direc_mat[i][j][4] = E

                elif max_score == t2:
                    # horizental
                    direc_mat[i][j][1]=1
                    direc_mat[i][j][3]=E
                else:
                    # vertical, seq1 has a gap
                    direc_mat[i][j][0]=1




def local_trace():
    # find the largest position(i,j)
    row_large =[]
    for i in range(len(score_mat)):
        row_max = max(score_mat[i])
        index = score_mat[i].index(row_max)
        row_large.append([row_max,index])

    t = [row_large[i][0] for i in range(len(row_large))]
    max_score = max(t)
    i = t.index(max_score)
    j = row_large[i][1]

    # i,j is the index of the largest element

    path_pair = []
    path_up = []
    path_down = []

    dest_i = i
    dest_j = j
    while(score_mat[i][j]!=0):
        score_mat[i][j] = 0
        tlist = direc_mat[i][j]
        if tlist[0]==1:
            path_up.append(seq1[i-1])
            path_down.append('-')
            i = i-1
        elif tlist[1]==1:
            path_down.append(seq2[j-1])
            path_up.append('-')
            j = j-1
        else:
            path_down.append(seq2[j-1])
            path_up.append(seq1[i-1])
            i = i -1
            j= j-1


def recal(i,j,mat=mat1):
    score_mat[i][j] = 0
    if i< len(seq1)+1:
        if direc_mat[i+1][j][0]==1:
            # do recalculate
            t1 = score_mat[i - 1][j - 1] + mat.get(seq1[i - 1], seq2[j - 1])

            # [s(i,j-n(gap)+g(gap]

            t2 = score_mat[i][j - 1] + direc_mat[i][j - 1][4]

            t3 = score_mat[i - 1][j] + direc_mat[i - 1][j][3]

            max_score = max(t1, t2, t3, 0)
            score_mat[i][j] = max_score

            if max_score > 0:
                if max_score == t1:
                    direc_mat[i][j][2] = 1
                    direc_mat[i][j][4] = E

                elif max_score == t2:
                    # horizental
                    direc_mat[i][j][1] = 1
                    direc_mat[i][j][3] = E
                else:
                    # vertical, seq1 has a gap
                    direc_mat[i][j][0] = 1
            recal(i+1,j)
    if j < len(seq2)+1:
        if direc_mat[i][j+1][1] ==1:
            t1 = score_mat[i - 1][j - 1] + mat.get(seq1[i - 1], seq2[j - 1])

            # [s(i,j-n(gap)+g(gap]

            t2 = score_mat[i][j - 1] + direc_mat[i][j - 1][4]

            t3 = score_mat[i - 1][j] + direc_mat[i - 1][j][3]

            max_score = max(t1, t2, t3, 0)
            score_mat[i][j] = max_score

            if max_score > 0:
                if max_score == t1:
                    direc_mat[i][j][2] = 1
                    direc_mat[i][j][4] = E

                elif max_score == t2:
                    # horizental
                    direc_mat[i][j][1] = 1
                    direc_mat[i][j][3] = E
                else:
                    # vertical, seq1 has a gap
                    direc_mat[i][j][0] = 1
                    recal(i,j+1)
    if i<len(seq1)+1 and len(seq2)+1:
        if direc_mat[i+1][j+1][2] ==1:










local_score(mat1)
print_mat2(score_mat)
local_trace()
print()
print_mat2(score_mat)



