from read_fasta import *
from Matrix import *
from numpy import *

base =-1.3

blosum = Matrix("BLOSUM62.txt")

amino = ['A', 'Q', 'L', 'S', 'R', 'E', 'K', 'T', 'N', 'G', 'M', 'W', 'D', 'H', 'F', 'Y', 'C', 'I', 'P', 'V','-']
# print(amino)

seq_aligned = read_fasta("WW-aligned-136.fasta")
seq2 =seq_aligned[0]
print(len(seq_aligned))

print(len(seq2))
for i in seq_aligned:
    pass
    #print(i)
fua = {}



def print_mat2(alist):
    for i in alist:
        for j in i:
            print('{0:4f}'.format(j), end=" ")
        print()


def fre(b, pos):
    count = 0
    for i in range(len(seq_aligned)):
        if b == seq_aligned[i][pos]:
            count += 1
    return count / len(seq_aligned)



def print_fua(dict):
    for key in dict:
        print(key, end=':')
        for j in dict[key]:
            print('{0:4f}'.format(j), end=" ")
        print()
    print()
    print()


# calculate the fua
for i in amino:
    t = []
    for j in range(len(seq_aligned[0])):
        t.append(fre(i, j))
    fua[i] = t
# print_fua(fua)
f = open("pa.txt")

#print_fua(fua)
# construct the pa dictionary
# what's the Pa for '-'?


pa = {}
for line in f:
    t = line.split()
    for i in range(0, len(t), 2):
        tt = t[i]
        pa[tt] = float(t[i + 1]) / 100

count =0
for i in seq_aligned:
    for j in i:
        if j=='-':
            count +=1
pa['-'] = count/(len(seq_aligned)*len(seq_aligned))
Nseq = len(seq_aligned)
alpa = Nseq - 1
beta = sqrt(Nseq)

# calculate the pua
pua = {}
for key in fua:
    t = []
    for i in range(len(seq_aligned[0])):
        t.append((alpa * fua[key][i] + beta * pa[key]) / (alpa + beta))

    pua[key] = t

# calculate the mua, the final PSSM matrix

mua = {}
for key in pua:
    t = []
    for i in range(len(seq_aligned[0])):
        t.append(log(pua[key][i] / pa[key]))
    mua[key] = t

gap1=[]
for i in range(len(seq_aligned[0])):
    gap1.append(fre('-',i))
print(gap1)
gap2 =[(tanh(i)-0.8) for i in gap1]
print(gap2)
print(pua['A'])
# print_fua(pua)
print_fua(mua)

f = open("mua.txt",'a')
for key in mua:
    s = ''+key+":"
    for j in mua[key]:
        s += '{0:4f}'.format(j)
    s += "\n"
    f.write(s)


# print the max value in the matrix
def max_value(dict=pua):
    result = []
    for i in range(len(seq_aligned[0])):
        t = [dict[key][i] for key in dict]
        result.append(amino[argmax(t)])
    return result


# print(max_value())






# begin score function


seq_list = read_fasta("protein-sequences.fasta")
seq1 = seq_list[0]
# initialize the score_mat and direc_mat
lseq2 = len(mua['A'])
score_mat = []
direc_mat = []

for i in range(len(seq1)+1):
    direc_mat.append([[0, 0, 0] for i in range(lseq2)])
for i in range(len(seq1)+1):
    score_mat.append([0 for i in range(lseq2)])


def print_mat2(alist):
    for i in alist:
        for j in i:
            print('{0:4f}'.format(j), end=" ")
        print()


def print_mat3(alist, n=3):
    for i in alist:
        for j in i:
            print(j[0:n], end='')
        print()


def local_pssm(istart=1, jstart=1):
    for i in range(istart, len(seq1)):
        for j in range(jstart,lseq2):
            t1 = score_mat[i - 1][j - 1] + mua[seq1[i]][j]

            t2 = score_mat[i][j - 1] + gap2[j-1]

            ## ??? doubt
            t3 = score_mat[i - 1][j] + gap2[j]

            max_score = max(t1, t2, t3, 0)
            score_mat[i][j] = max_score

            if max_score > 0:
                if max_score == t1:
                    direc_mat[i][j][2] = 1
                if max_score == t2:
                    direc_mat[i][j][1] = 1
                if max_score == t3:
                    direc_mat[i][j][0] = 1

def print_pathpair(alist):
    for path in alist:
        mark1=[]
        mark2=[]
        for i in range(len(path[0])):
            if path[0][i]!='-' and path[1][i]!='-':
                if path[0][i] == path[1][i]:
                    mark1.append('.')
                    mark2.append('.')
                elif mua[path[0][i]][i]>0:
                    mark1.append('.')
                    mark2.append(' ')
                else:
                    mark1.append(' ')
                    mark2.append(' ')
            else:
                mark1.append(' ')
                mark2.append(' ')

        path_u = ''.join(path[0])
        path_m = ''.join(mark1)
        path_2 = ''.join(mark2)

        path_d = ''.join(path[1])


        print(path_u)
        print(path_m)
        print(path_2)
        print(path_d)



def traceback(k=1):
    path_pair = []
    path_down = []
    path_up=[]

    # find the largest element in the j column
    t = [score_mat[i][lseq2-1] for i in range(len(seq1))]
    i = argmax(t)

    recal_pair=[]


    print('score:',score_mat[i][lseq2-1])
    print("the coordinate of the starting tracing coordinate:", i)

    queue = []

    queue.append([path_up[0:len(path_up)], path_down[0:len(path_down)], i, lseq2-1])
    while (len(queue) > 0):

        t = queue.pop(0)
        i = t[2]
        j = t[3]
        path_up = t[0]
        path_down = t[1]


        if score_mat[i][j]==0:
            path_pair.append([path_up,path_down])
        # scan all possible path and append it to the queue
        else:
            if direc_mat[i][j][1] == 1:
                path_up.insert(0, '-')
                j = j - 1

                path_down.insert(0, '*')

                recal_pair.append([i,j])

                if len(queue) < k:
                    queue.append([path_up[0:len(path_up)], path_down[0:len(path_down)], i, j])
                j = j + 1
                path_up = path_up[1:]
                path_down = path_down[1:]
            if direc_mat[i][j][0] == 1:
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
                path_down.insert(0, '*')

                recal_pair.append([i,j])

                if len(queue) < k:
                    queue.append([path_up[0:len(path_up)], path_down[0:len(path_down)], i, j])
                i = i + 1
                j = j + 1
                path_up = path_up[1:]
                path_down = path_down[1:]

            score_mat[i][j]=0

    #print(path_pair)

    print("recal_pair:",recal_pair)

    for k in range(len(recal_pair) - 1, -1, -1):
        t = recal_pair[k]
        recal(recal_pair[k][0], recal_pair[k][1])

    return path_pair




def print_pathpair(alist):
    for path in alist:
        path_u = ''.join(path[0])
        path_d = ''.join(path[1])

        print(path_u)
        print(path_d)


def recal(i,j):
    score_mat[i][j] = 0
    if i< len(seq1)-2 and j <lseq2-2:
        if direc_mat[i+1][j][0]==0 and direc_mat[i][j+1][1]==0 and direc_mat[i+1][j+1][2]==0:
            return
    if i< len(seq1)-2:
        i = i+1
        if direc_mat[i][j][0]==1:
            # do recalculate
            t1 = score_mat[i - 1][j - 1] + mua[seq1[i]][j]

            t2 = score_mat[i][j - 1] + gap2[j - 1]

            ## ??? doubt
            t3 = score_mat[i - 1][j] + gap2[j]
            max_score = max(t1, t2, t3, 0)
            score_mat[i][j] = max_score

            if max_score > 0:
                if max_score == t1:
                    direc_mat[i][j][2] = 1


                elif max_score == t2:
                    # horizental
                    direc_mat[i][j][1] = 1

                else:
                    # vertical, seq1 has a gap
                    direc_mat[i][j][0] = 1
            else:
                direc_mat[i][j][0]=0
                direc_mat[i][j][1]=0
                direc_mat[i][j][2]=0
            recal(i,j)


        i = i-1

    if j < lseq2-2:
        j = j+1
        if direc_mat[i][j][1] ==1:
            t1 = score_mat[i - 1][j - 1] + mua[seq1[i]][j]

            t2 = score_mat[i][j - 1] + gap2[j - 1]

            ## ??? doubt
            t3 = score_mat[i - 1][j] + gap2[j]

            max_score = max(t1, t2, t3, 0)
            score_mat[i][j] = max_score

            if max_score > 0:
                if max_score == t1:
                    direc_mat[i][j][2] = 1


                elif max_score == t2:
                    # horizental
                    direc_mat[i][j][1] = 1

                else:
                    # vertical, seq1 has a gap
                    direc_mat[i][j][0] = 1
            else:
                direc_mat[i][j][0]=0
                direc_mat[i][j][1]=0
                direc_mat[i][j][2]=0
            recal(i, j)
        j = j-1

    if i<len(seq1)-2 and j<lseq2-1:
        i = i+1
        j = j+1
        if direc_mat[i][j][2] ==1:
            t1 = score_mat[i - 1][j - 1] + mua[seq1[i]][j]

            t2 = score_mat[i][j - 1] + gap2[j - 1]

            ## ??? doubt
            t3 = score_mat[i - 1][j] + gap2[j]

            max_score = max(t1, t2, t3, 0)
            score_mat[i][j] = max_score

            if max_score > 0:
                if max_score == t1:
                    direc_mat[i][j][2] = 1

                elif max_score == t2:
                    # horizental
                    direc_mat[i][j][1] = 1
                else:
                    # vertical, seq1 has a gap
                    direc_mat[i][j][0] = 1
            else:
                direc_mat[i][j][0]=0
                direc_mat[i][j][1]=0
                direc_mat[i][j][2]=0
            recal(i, j)


        i = i -1
        j = j -1






local_pssm()
print_mat3(direc_mat,3)

#print(lseq2)
#print(len(seq1))
# print_mat3(direc_mat)
print_pathpair(traceback(2))
print_pathpair(traceback(2))







