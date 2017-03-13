from read_fasta import *
from Matrix import *
from numpy import *
blosum = Matrix("BLOSUM62.txt")



amino = ['A', 'Q', 'L', 'S', 'R', 'E', 'K', 'T', 'N', 'G', 'M', 'W', 'D', 'H', 'F', 'Y', 'C', 'I', 'P', 'V','-']
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
    for key in dict:
        print(key,end=':')
        for j in dict[key]:
            print('{0:4f}'.format(j), end=" ")
        print()
    print()
    print()


for i in amino:
    t =[]
    for j in range(len(seq_aligned[0])):
        t.append(fre(i,j))
    fua[i]=t
print_fua(fua)
f = open("pa.txt")

pa = {}
pa['-']=1
for line in f:
    t = line.split()
    for i in range(0,len(t),2):
        tt=t[i]
        pa[tt] = float(t[i+1])/100


Nseq = len(seq_aligned)
alpa = Nseq-1
beta = sqrt(Nseq)

print(pa.keys())
pua = {}
for key in fua:
    t =[]
    for i in range(len(seq_aligned[0])):
       t.append((alpa*fua[key][i]+beta*pa[key])/(alpa+beta))

    pua[key] =t
mua = {}
for key in pua:
    t =[]
    for i in range(len(seq_aligned[0])):
        t.append(log(pua[key][i]/pa[key]))
    mua[key] = t
print_fua(pua)
print_fua(mua)


# print the max value in the matrix
def max_value(dict=pua,num =1):
    result = []
    for i in range(len(seq_aligned[0])):
        t = [dict[key][i] for key in dict]
        result.append(amino[argmax(t)])
    return result

print(max_value())

