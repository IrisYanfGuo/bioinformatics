from read_fasta import *
from Matrix import *
from numpy import *

blosum = Matrix("BLOSUM62.txt")

amino = ['A', 'Q', 'L', 'S', 'R', 'E', 'K', 'T', 'N', 'G', 'M', 'W', 'D', 'H', 'F', 'Y', 'C', 'I', 'P', 'V', '-']
# print(amino)

seq_aligned = read_fasta("WW-aligned-136.fasta")
for i in seq_aligned:
    print(i)
fua = {}


def print_mat2(alist):
    for i in alist:
        for j in i:
            print('{0:4d}'.format(j), end=" ")
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

# construct the pa dictionary
# what's the Pa for '-'?
pa = {}
pa['-'] = 1
for line in f:
    t = line.split()
    for i in range(0, len(t), 2):
        tt = t[i]
        pa[tt] = float(t[i + 1]) / 100

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


# print_fua(pua)
print_fua(mua)


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
lseq2 = len(mua['A'])
# initialize the score_mat and direc_mat
score_mat = []
direc_mat = []

for i in range(len(seq11) + 1):
    direc_mat.append([[0, 0, 0] for i in range(lseq2 + 1)])
for i in range(len(seq11) + 1):
    score_mat.append([0 for i in range(lseq2 + 1)])


def print_mat2(alist):
    for i in alist:
        for j in i:
            print('{0:4d}'.format(j), end=" ")
        print()


def print_mat3(alist, n=3):
    for i in alist:
        for j in i:
            print(j[0:n], end='')
        print()


def local_pssm(istart=1, jstart=1):
    for i in range(istart, len(seq1) + 1):
        for j in range(jstart, lseq2 + 1):
            t1 = score_mat[i - 1][j - 1] + mua[seq1[i]][j]

            t2 = score_mat[i][j - 1] + mua['-'][j - 1]

            t3 = score_mat[i - 1][j] + mua['-'][j]

            max_score = max(t1, t2, t3, 0)
            score_mat[i][j] = max_score

            if max_score > 0:
                if max_score == t1:
                    direc_mat[i][j][2] = 1
                if max_score == t2:
                    direc_mat[i][j][1] = 1
                if max_score == t3:
                    direc_mat[i][j][0] = 1


def local_trace(k=5):
    # find the largest position(i,j)
    row_large = []
    for i in range(len(score_mat)):
        row_max = max(score_mat[i])
        index = score_mat[i].index(row_max)
        row_large.append([row_max, index])

    t = [row_large[i][0] for i in range(len(row_large))]
    max_score = max(t)
    i = t.index(max_score)
    j = row_large[i][1]

    path = []
    path_seq2 = []
    recal_pair=[]

    queue = []
    queue.append([path_seq2[0:len([path_seq2])],i,j])

    print("score=", score_mat[i][j])

    while(len(queue)>0):
        t = queue.pop(0)
        path_seq2 = t[0]
        i = t[1]
        j = t[2]

        if score_mat[i][j] ==0:
            path.append(path_seq2)
        else:
            if direc_mat[i][j][0]==1:
                recal_pair.append(i,j)
                if len(queue)<k :
                    queue.append([path_seq2[0:len(path_seq2)],i,j])
            if direc_mat[i][j][1]==1 :
                recal_pair.append(i,j)

                path_seq2








