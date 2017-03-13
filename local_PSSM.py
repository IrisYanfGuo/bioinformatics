# in a non recursive way
from Matrix2 import *
from PSSM import *
# affine gap penalty
E = -4
I = -4 # affine gap


seq_list = read_fasta("protein-sequences.fasta")


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




def local_trace(k=5):
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
                j = j + 1
            if direc_mat[i][j][1] == 1:
                path_down.insert(0, '-')
                i = i - 1
                path_up.insert(0, seq1[i])
                recal_pair.append([i,j])
                if len(queue) < k:
                    queue.append([path_up[0:len(path_up)], path_down[0:len(path_down)], i, j])
                i = i + 1

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

        for i in range(len(path_pair)):
            print("path ",i+1,":")
            print(''.join(path_pair[i][0]))
            print(''.join(path_pair[i][1]))
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
local_trace()

#print_mat2(score_mat)





