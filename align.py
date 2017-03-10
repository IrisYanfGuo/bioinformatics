# real things
from Matrix import *
from read_fasta import *

seq_list = read_fasta("WW-sequence.fasta")
sequ1 = seq_list[0][0:6]
sequ2 = seq_list[1][0:5]
mat1 = Matrix("BLOSUM50.txt")
print(sequ1, sequ2)


# mat1.print()



def memorize(f):
    result = None

    def wrapper(a1, a2, a3):
        nonlocal result
        if result is None:
            result = f(a1, a2, a3)

    return wrapper

    # initialize the score mat


score_mat = []

for i in range(len(sequ1) + 1):
    score_mat.append(['#' for i in range(len(sequ2) + 1)])
score_mat[0][0] = 0
t = 0
for i in range(1, len(sequ2) + 1):
    t += mat1.get('*', sequ2[i - 1])
    score_mat[0][i] = t
t = 0
for i in range(1, len(sequ1) + 1):
    t += mat1.get('*', sequ1[i - 1])
    score_mat[i][0] = t
# print(score_mat)

# create the direction mat
direc_mat = []
for i in range(len(sequ1) + 1):
    direc_mat.append([[0, 0, 0] for i in range(len(sequ2) + 1)])
# [0,1,0] means horizontal
# [1,0,0] means vertical
# [0,0,1] means diagonal
for i in range(1, len(sequ2) + 1):
    direc_mat[0][i] = [0, 1, 0]
for i in range(1, len(sequ1) + 1):
    direc_mat[i][0] = [1, 0, 0]
    # print(direc_mat)


def score(seq1, seq2, mat, pos1, pos2):

    if (score_mat[pos1][pos2] != '#'):
        #print(score_mat[pos1][pos2],pos1,pos2)
        return score_mat[pos1][pos2]
    else:
        # insert a gap on seq2
        t1 = score(seq1, seq2, mat, pos1 - 1, pos2) + mat.get(seq1[pos1 - 1], '*')
        # insert a gap on seq1
        t2 = score(seq1, seq2, mat, pos1, pos2 - 1) + mat.get(seq2[pos2 - 1], '*')

        t3 = score(seq1, seq2, mat, pos1 - 1, pos2 - 1) + mat.get(seq1[pos1 - 1], seq2[pos2 - 1])

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


print(score(sequ1, sequ2, mat1, len(sequ1), len(sequ2)))

def print_mat2(alist):
    for i in alist:
        for j in i:
            print('{0:4d}'.format(j),end=" ")
        print()

def print_mat3(alist):
    for i in alist:
        print(i)


def traceback(dir_m,k,seq1,seq2):
    allPath = []
    seq1_len = len(seq1)
    seq2_len = len(seq2)



print_mat3(direc_mat)
print_mat2(score_mat)