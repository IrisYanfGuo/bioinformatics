from pre_process import *
from numpy import log
import matplotlib.pyplot as plt

dssp2 = pre_process("dssp_info.txt", "dssp.txt")
stride = pre_process("stride_info.txt", "stride.txt")

fsr, fs, fr = fSR(dssp2)
# print(fs)
# print(fr)
total = fs['C'] + fs['E'] + fs['C']

# do some pre_process, read the protein sequence and structure in to a 2d list
f = open("dssp_protein")

list = f.readlines()
list2 = []
for i in range(0, len(list), 3):
    temp = []
    temp.append(list[i + 1].strip())
    temp.append(list[i + 2].strip())
    list2.append(temp)
f.close()
print(list2)
# dealing the statistics
fsrm = {}
# initialize
amino = []
for i in amino_dict.keys():
    amino.append(amino_dict[i])
amino.append('?')
print(amino)
for i in ['C', 'E', 'H']:
    temp3 = {}
    for j in amino:
        temp2 = {}
        for k in amino:
            temp = {}
            for n in range(1, 9):
                temp[n] = 1

            temp2[k] = temp

        temp3[j] = temp2
    fsrm[i] = temp3

for i in list2:
    for j in range(len(i[0])):
        for k in range(1, 9):

            fsrm[i[1][j]][i[0][j]][i[0][k]][k] += 1

for i in fsrm.keys():
    print(fsrm[i])


def gor3(alist):
    result = []
    for i in range(len(alist)):

        helix = log(fsr['H'][alist[i]] / (fsr['C'][alist[i]] + fsr['E'][alist[i]])) + log(
            (fs['E'] + fs['C']) / fs['H'])
        coil = log(fsr['C'][alist[i]] / (fsr['H'][alist[i]] + fsr['E'][alist[i]])) + log(
            (fs['H'] + fs['E']) / fs['C'])
        sheet = log(fsr['E'][alist[i]] / (fsr['C'][alist[i]] + fsr['H'][alist[i]])) + log(
            (fs['H'] + fs['C']) / fs['E'])


        for j in range(-8, 9):
            t = i + j
            if t > 0 and t < len(alist):
                if j != 0:
                    helix += log(fsrm['H'][alist[i]][alist[t]][abs(j)] / (
                    fsrm['E'][alist[i]][alist[t]][abs(j)] + fsrm['C'][alist[i]][alist[t]][abs(j)])) + log(
                        (fsr['E'][alist[i]] + fsr['C'][alist[i]]) / fsr['H'][alist[i]])

                    coil += log(fsrm['C'][alist[i]][alist[t]][abs(j)] / (
                        fsrm['E'][alist[i]][alist[t]][abs(j)] + fsrm['H'][alist[i]][alist[t]][abs(j)])) + log(
                        (fsr['E'][alist[i]] + fsr['H'][alist[i]]) / fsr['C'][alist[i]])
                    sheet += log(fsrm['E'][alist[i]][alist[t]][abs(j)] / (
                        fsrm['H'][alist[i]][alist[t]][abs(j)] + fsrm['C'][alist[i]][alist[t]][abs(j)])) + log(
                        (fsr['H'][alist[i]] + fsr['C'][alist[i]]) / fsr['E'][alist[i]])

        if coil == max(sheet, coil, helix):
            result.append('C')
        elif sheet == max(sheet, coil, helix):
            result.append('E')

        else:
            result.append('H')



    return ''.join(result)



'''
f = open('dssp_protein')
protein=[]
stru =[]

t = f.readlines()
for i in range(0,len(t),3):
    protein.append(t[i+1].strip())
    stru.append(t[i+2].strip())
f.close()

f = open("aa.txt",'w')
for i in protein:
    f.write(gor3(i))
f.write('\n')
for i in stru:
    f.write(i)
f.close()
'''

'''
f = open('dssp_protein')
protein=[]
stru =[]

t = f.readlines()
for i in range(0,len(t),3):
    protein.append(t[i+1].strip())
    stru.append(t[i+2].strip())
f.close()


aaa = 'G?NQSIIFTEQLTWDVQLSAIHFTAQQQG?VIDCYIGQKVLEHLAAEKINNSEQALSLFEQFRFDIEEQAEKLIEQEAFDVQGHIQVERVD'
bbb ='CCCCCEEECCCCEEECCCCEEEEEEEECCEEEEEEEEHHHHHHHHCCCCCCHHHHHHHHHHCHHHHHHHHHHHHHCCCCCCCCCEEECCCC'
f = open("aa.txt",'w')
f.write(gor3(aaa))
f.write('\n')
f.write(bbb)

'''