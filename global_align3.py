# in a non recursive way
from Matrix import *
from read_fasta import *

# affine gap penalty
E = -4
I = -5 # affine gap

seq_list = read_fasta("protein-sequences.fasta")
seq_list2 = read_fasta("WW-sequence.fasta")
seq111 =seq_list2[0]
print(len(seq111))

seq2222 =seq_list2[2]
print(len(seq2222))
sequ1 = 'THISLINE'
sequ2 = 'ISALIGNED'
mat1 = Matrix("BLOSUM62.txt")



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


def traceback(direc_mat,seq1,seq2,k=4):
    path_pair = []
    path_up = []
    path_down = []
    mark1 =[]
    mark2 =[]

    # use queue to trace multiple path

    queue = []
    queue.append([path_up[0:len(path_up)], path_down[0:len(path_down)], len(seq1), len(seq2)])
    while(len(queue)>0):
        t = queue.pop(0)
        i = t[2]
        j = t[3]
        path_up = t[0]
        path_down=t[1]

        if i ==0 or j ==0:
            if j==0 and i ==0 :
                path_pair.append([path_up,path_down])
            elif i ==0:
                for m in range(j+1,-1,-1):
                    path_up.insert(0,'-')
                    path_down.insert(0,seq2[m])
                path_pair.append([path_up,path_down])
            else:
                for m in range(i+1,-1,-1):
                    path_down.insert(0,'-')
                    path_up.insert(0,seq1[m])
                path_pair.append([path_up,path_down])
        #scan all possible path and append it to the queue
        else:
            if direc_mat[i][j][0] == 1:
                path_up.insert(0, '-')
                j = j - 1

                path_down.insert(0, seq2[j])
                if len(queue)<k :
                    queue.append([path_up[0:len(path_up)], path_down[0:len(path_down)], i, j])
                j = j+1
                path_up = path_up[1:]
                path_down = path_down[1:]
            if direc_mat[i][j][1] == 1:
                path_down.insert(0, '-')
                i = i - 1
                path_up.insert(0, seq1[i])
                if len(queue)<k :
                    queue.append([path_up[0:len(path_up)], path_down[0:len(path_down)], i, j])
                i = i+1
                path_up = path_up[1:]
                path_down = path_down[1:]

            if direc_mat[i][j][2] == 1:
                i = i - 1
                j = j - 1
                path_up.insert(0, seq1[i])
                path_down.insert(0, seq2[j])

                if len(queue)<k :
                    queue.append([path_up[0:len(path_up)], path_down[0:len(path_down)], i, j])
                i = i+1
                j = j+1
                path_up = path_up[1:]
                path_down = path_down[1:]



    return path_pair


def print_pathpair(alist):
    for path in alist:
        mark1=[]
        mark2=[]
        for i in range(len(path[0])):
            if path[0][i]!='-' and path[1][i]!='-':
                if path[0][i] == path[1][i]:
                    mark1.append(':')

                elif mat1.get(path[0][i],path[1][i])>0:
                    mark1.append('.')

                else:
                    mark1.append(' ')

            else:
                mark1.append(' ')


        path_u = ''.join(path[0])
        path_m = ''.join(mark1)


        path_d = ''.join(path[1])


        print(path_u)
        print(path_m)
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

score2()

fun(seq111,seq2222)
print(seq2222)