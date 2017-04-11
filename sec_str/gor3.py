from pre_process import *
from numpy import log
import matplotlib.pyplot as plt

dssp2 = pre_process("dssp_info.txt", "dssp.txt")
stride = pre_process("stride_info.txt", "stride.txt")

fsr,fs,fr= fSR(dssp2)
#print(fs)
#print(fr)
total = fs['C']+fs['E']+fs['C']



# do some pre_process, read the protein sequence and structure in to a 2d list
f = open("dssp_protein")

list = f.readlines()
list2 =[]
for i in range(0,len(list),3):
    temp =[]
    temp.append(list[i+1].strip())
    temp.append(list[i+2].strip())
    list2.append(temp)
f.close()
# dealing the statistics
fsrm={}
# initialize
amino =[]
for i in amino_dict.keys():
    amino.append(amino_dict[i])
amino.append('?')
print(amino)
for i in ['C','E','H']:
    temp3 = {}
    for j in amino:
        temp2 = {}
        for k in amino:
            temp = {}
            for n in range(1,9):

                temp[n]=0

            temp2[k] = temp

        temp3[j] = temp2
    fsrm[i] = temp3

for i in list2:
    for j in range(len(i[0])):
        for k in range(1,9):
            fsrm[i[1][j]][i[0][j]][i[0][k]][k] +=1

for i in fsrm.keys():
    print(fsrm[i])
