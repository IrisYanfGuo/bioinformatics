from read_fasta import *
from Matrix import *
blosum = Matrix("BLOSUM62.txt")
blosum.print()

amino = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '-']
print(amino)




seq_aligned = read_fasta("WW-aligned-136.fasta")


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
    return count/len(seq_aligned)

def print_fua(dict):
    for key in fua:
        print(key,end=':')
        for j in fua[key]:
            print('{0:4f}'.format(j), end=" ")
        print()


for i in amino:
    t =[]
    for j in range(len(seq_aligned)):
        t.append(fre(i,j))
    fua[i]=t


print_fua(fua)



