# in a non recursive way
from Matrix import *
from read_fasta import *

# affine gap penalty
E = -4
I = -6 # affine gap

seq_list = read_fasta("WW-sequence.fasta")
seq1 =seq_list[0][1:30]
seq2 =seq_list[1][1:20]
sequ1 = 'THISLINE'
sequ2 = 'ISALIGNED'
slen1 = len(seq1)
slen2 = len(seq2)
mat1 = Matrix("BLOSUM50.txt")
print(seq1, seq2)

score_mat = []
direc_mat = []
for i in range(len(seq1) + 1):
    # default, opening gap
    direc_mat.append([[0, 0, 0,I,I] for i in range(len(seq2) + 1)])

for i in range(len(seq1) + 1):
    score_mat.append([0 for i in range(len(seq2) + 1)])
# initialize
score_mat[0][0] = 0

direc_mat[0][1] = [0,1,0,I,E]
score_mat[0][1]=I
direc_mat[1][0]=[1,0,0,E,I]
score_mat[1][0]=I

for i in range(2, len(seq2)+1):
    direc_mat[0][i] = [0,1,0,I,E]
    score_mat[0][i] = score_mat[0][i-1] + direc_mat[0][i-1][4]
for i in range(2, len(seq1)+1):
    direc_mat[i][0] = [1,0,0,E,I]
    score_mat[i][0] = score_mat[i-1][0]+ direc_mat[i-1][0][3]
# create the direction mat


# [0,1,0] means seq 2 has a gap
# [1,0,0] means seq1 has a gap
# [0,0,1] means diagonal

    # print(direc_mat)


def score2(mat):
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):

            t1 = score_mat[i - 1][j] + direc_mat[i-1][j][3]

            t2 = score_mat[i][j - 1] +direc_mat[i][j-1][4]

            t3 = score_mat[i - 1][j - 1] + mat.get(seq1[i - 1], seq2[j - 1])

            maxscore = max(t1, t2, t3)
            score_mat[i][j] = maxscore
            if (t1 == maxscore):
                direc_mat[i][j][1] = 1
                # affine gap
                direc_mat[i][j][3] = E
            if (t2 == maxscore):
                direc_mat[i][j][0] = 1
                direc_mat[i][j][4] = E
            if (t3 == maxscore):
                direc_mat[i][j][2] = 1
    return score_mat[len(seq1)][len(seq2)]


def print_mat2(alist):
    for i in alist:
        for j in i:
            print('{0:4d}'.format(j), end=" ")
        print()


def print_mat3(alist):
    for i in alist:
        print(i)


def traceback():
    i = slen1
    j = slen2
    path_pair = []
    path_up = []
    path_down = []

    # use queue to trace multiple path
    queue = []
    queue.append([slen1, slen2])
    while len(queue) != 0:
        t = queue.pop(0)
        i = t[0]
        j = t[1]
        while (i > 0 and j > 0):
            tlist = direc_mat[i][j]
            if tlist[0] == 1:
                path_up.append('-')
                j = j - 1
                path_down.append(seq2[j])
            elif tlist[1] == 1:
                path_down.append('-')
                i = i - 1
                path_up.append(seq1[i])
            elif tlist[2] == 1:
                i = i - 1
                j = j - 1
                path_up.append(seq1[i])
                path_down.append(seq2[j])
        path_pair.append([path_up, path_down])
    return path_pair


def print_pathpair(alist):
    for path in alist:
        path_u = ''.join(path[0])
        path_d = ''.join(path[1])
        print(path_u)
        print(path_d)


print(score2(mat1))
print_pathpair(traceback())
#print_mat2(score_mat)
#print_mat3(direc_mat)
print_mat2(score_mat)
print_mat3(direc_mat)



