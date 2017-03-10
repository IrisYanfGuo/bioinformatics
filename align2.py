from Matrix import *
from read_fasta import *

seq_list = read_fasta("WW-sequence.fasta")
seq1 = seq_list[0][:20]
seq2 = seq_list[1][:15]
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


def score(mat, pos1, pos2):
    if (score_mat[pos1][pos2] != '#'):
        # print(score_mat[pos1][pos2],pos1,pos2)
        return score_mat[pos1][pos2]
    else:
        # insert a gap on seq2
        t1 = score(mat, pos1 - 1, pos2) + mat.get(seq1[pos1 - 1], '*')
        # insert a gap on seq1
        t2 = score(mat, pos1, pos2 - 1) + mat.get(seq2[pos2 - 1], '*')

        t3 = score(mat, pos1 - 1, pos2 - 1) + mat.get(seq1[pos1 - 1], seq2[pos2 - 1])

        # compare them
        maxscore = max(t1, t2, t3)
        score_mat[pos1][pos2] = maxscore
        if (t1 == maxscore):
            direc_mat[pos1][pos2][1] = 1
        if (t2 == maxscore):
            direc_mat[pos1][pos2][0] = 1
        if (t3 == maxscore):
            direc_mat[pos1][pos2][2] = 1
        return maxscore


print(score(mat1, len(seq1), len(seq2)))


def print_mat2(alist):
    for i in alist:
        for j in i:
            print('{0:4d}'.format(j), end=" ")
        print()


def print_mat3(alist):
    for i in alist:
        print(i)


print_mat2(score_mat)
print_mat3(direc_mat)


# [0,1,0] means seq 2 has a gap
# [1,0,0] means seq1 has a gap
# [0,0,1] means diagonal

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
        while (i != 0 or j != 0):
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


print_pathpair(traceback())
