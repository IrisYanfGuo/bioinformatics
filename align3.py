# in a non recursive way
from Matrix import *
from read_fasta import *

# affine gap penalty
E = -4
I = -6 # affine gap

seq_list = read_fasta("WW-homo-136.fasta")
seq_list2 = read_fasta("WW-sequence.fasta")
seq111 =seq_list[0]
seq2222 =seq_list[1][:-4]
sequ1 = 'THISLINE'
sequ2 = 'ISALIGNED'
mat1 = Matrix("BLOSUM50.txt")



def score2(mat,seq1,seq2,score_mat,direc_mat):
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


def traceback(direc_mat,seq1,seq2):
    path_pair = []
    path_up = []
    path_down = []
    mark1 =[]
    mark2 =[]

    # use queue to trace multiple path




    i = len(seq1)
    j = len(seq2)


    while (i > 0 and j > 0):
        tlist = direc_mat[i][j]
        if tlist[0] == 1:
            path_up.insert(0, '-')
            mark1.insert(0,' ')
            mark2.insert(0,' ')
            j = j - 1
            path_down.insert(0, seq2[j])
        elif tlist[1] == 1:
            path_down.insert(0, '-')
            mark1.insert(0, ' ')
            mark2.insert(0, ' ')
            i = i - 1
            path_up.insert(0, seq1[i])
        elif tlist[2] == 1:
            i = i - 1
            j = j - 1
            if seq1[i]==seq2[j]:
                mark1.insert(0,'.')
                mark2.insert(0,'.')
            elif mat1.get(seq1[i],seq2[j])>0:
                mark1.insert(0,' ')
                mark2.insert(0,'.')
            else:
                mark1.insert(0,' ')
                mark2.insert(0,' ')
            path_up.insert(0, seq1[i])
            path_down.insert(0, seq2[j])
    path_pair.append([path_up, mark1,mark2,path_down])
    return path_pair


def print_pathpair(alist):
    for path in alist:
        path_u = ''.join(path[0])
        path_d = ''.join(path[3])
        mark1=''.join(path[1])
        mark2 = ''.join(path[2])

        print(path_u)
        print(mark1)
        print(mark2)
        print(path_d)


#print_mat3(direc_mat)
#print_mat2(score_mat)
#print_mat3(direc_mat)

def fun(seq1,seq2):
    score_mat = []
    direc_mat = []
    for i in range(len(seq1) + 1):
        # default, opening gap
        direc_mat.append([[0, 0, 0, I, I] for i in range(len(seq2) + 1)])

    for i in range(len(seq1) + 1):
        score_mat.append([0 for i in range(len(seq2) + 1)])
    # initialize
    score_mat[0][0] = 0

    direc_mat[0][1] = [0, 1, 0, I, E]
    score_mat[0][1] = I
    direc_mat[1][0] = [1, 0, 0, E, I]
    score_mat[1][0] = I

    for i in range(2, len(seq2) + 1):
        direc_mat[0][i] = [0, 1, 0, I, E]
        score_mat[0][i] = score_mat[0][i - 1] + direc_mat[0][i - 1][4]
    for i in range(2, len(seq1) + 1):
        direc_mat[i][0] = [1, 0, 0, E, I]
        score_mat[i][0] = score_mat[i - 1][0] + direc_mat[i - 1][0][3]

    print("score:",score2(mat1,seq1,seq2,score_mat,direc_mat))
    print_pathpair(traceback(direc_mat,seq1,seq2))
    print()





for i in seq_list:
    for j in seq_list:
        if i!=j:
            fun(i,j)