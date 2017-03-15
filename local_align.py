# in a non recursive way
from Matrix import *
from read_fasta import *

# affine gap penalty
E = -4
I = -6 # affine gap


seq_list = read_fasta("WW-homo-136.fasta")
seq_list2 = read_fasta("protein-sequences.fasta")

seq1 =seq_list[3]
seq2= seq_list2[1]
seq11= 'THISLINE'
seq22= 'ISALIGNED'
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


def print_mat3(alist,n=3):
    for i in alist:
        for j in i:
            print(j[0:n],end='')
        print()


def local_score(mat):
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):

            t1 = score_mat[i - 1][j - 1] + mat.get(seq1[i - 1], seq2[j - 1])

            # [s(i,j-n(gap)+g(gap]

            t2 =score_mat[i][j-1] +direc_mat[i][j-1][4]



            t3 = score_mat[i-1][j] + direc_mat[i-1][j][3]



            max_score = max(t1, t2, t3, 0)
            score_mat[i][j] = max_score

            if max_score > 0:
                if max_score == t1:
                    direc_mat[i][j][2] = 1


                elif max_score == t2:
                    # horizental
                    direc_mat[i][j][1]=1
                    direc_mat[i][j][4]=E
                else:
                    # vertical, seq1 has a gap
                    direc_mat[i][j][0]=1
                    direc_mat[i][j][3] = E


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





def local_trace(k=3):
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

    print(i,j)

    # i,j is the index of the largest element

    path_pair = []
    path_up = []
    path_down = []

    #path_up.insert(0,str(i))
    #path_down.insert(0,str(j))

    queue = []
    queue.append([path_up[0:len(path_up)], path_down[0:len(path_down)], i, j])

    recal_pair=[]

    print("score=",score_mat[i][j])
    while (len(queue) > 0):
        t = queue.pop(0)
        i = t[2]
        j = t[3]
        path_up = t[0]
        path_down = t[1]

        if score_mat[i][j] == 0:
            path_pair.append([path_up,path_down])
        # scan all possible path and append it to the queue
        else:
            if direc_mat[i][j][0] == 1:
                path_up.insert(0, '-')
                j = j - 1

                path_down.insert(0, seq2[j])
                recal_pair.append([i,j])
                if len(queue) < k:
                    queue.append([path_up[0:len(path_up)], path_down[0:len(path_down)], i, j])
                path_up = path_up[1:]
                path_down = path_down[1:]
                j = j + 1
            if direc_mat[i][j][1] == 1:
                path_down.insert(0, '-')
                i = i - 1
                path_up.insert(0, seq1[i])
                recal_pair.append([i,j])
                if len(queue) < k:
                    queue.append([path_up[0:len(path_up)], path_down[0:len(path_down)], i, j])
                i = i + 1
                path_up = path_up[1:]
                path_down = path_down[1:]

            if direc_mat[i][j][2] == 1:
                i = i - 1
                j = j - 1
                path_up.insert(0, seq1[i])
                path_down.insert(0, seq2[j])
                recal_pair.append([i,j])
                if len(queue) < k:
                    queue.append([path_up[0:len(path_up)], path_down[0:len(path_down)], i, j])
                i = i + 1
                j = j + 1
                path_up = path_up[1:]
                path_down = path_down[1:]

        print_pathpair(path_pair)

   # print("recal_pair:",recal_pair)

    for k in range(len(recal_pair)-1,-1,-1):
        t=recal_pair[k]
        recal(recal_pair[k][0],recal_pair[k][1])



def recal(i,j,mat=mat1):
    score_mat[i][j] = 0
    if i< len(seq1) and j <len(seq2):
        if direc_mat[i+1][j][0]==0 and direc_mat[i][j+1][1]==0 and direc_mat[i+1][j+1][2]==0:
            return
    if i< len(seq1):
        i = i+1
        if direc_mat[i][j][0]==1:
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
            recal(i,j)


        i = i-1

    if j < len(seq2):
        j = j+1
        if direc_mat[i][j][1] ==1:
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
            recal(i, j)
        j = j-1

    if i<len(seq1) and j<len(seq2):
        i = i+1
        j = j+1
        if direc_mat[i][j][2] ==1:
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
            recal(i, j)


        i = i -1
        j = j -1











local_score(mat1)
#print_mat2(score_mat)
#print_mat3(direc_mat)
local_trace()
local_trace(4)
local_trace()

#print_mat2(score_mat)





