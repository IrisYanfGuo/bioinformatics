from pre_process import *
from numpy import log
import time
from numpy import sqrt
import matplotlib.pyplot as plt
'''
t1 = time.time()
dssp2 = pre_process("./txtfile/dssp_info.txt", "./txtfile/dssp.txt")
stride = pre_process("./txtfile/stride_info.txt", "./txtfile/stride.txt")

fsr, fs, fr = fSR(dssp2)
# print(fs)
# print(fr)

total = fs['C'] + fs['E'] + fs['C']

# do some pre_process, read the protein sequence and structure in to a 2d list
f = open("./txtfile/dssp_protein")

list = f.readlines()
list2 = []
for i in range(0, len(list), 3):
    temp = []
    temp.append(list[i + 1].strip())
    temp.append(list[i + 2].strip())
    temp.append(list[i].strip())
    list2.append(temp)
f.close()
# print(list2)
# dealing the statistics
fsrm = {}
# initialize
amino = []
for i in amino_dict.keys():
    amino.append(amino_dict[i])
amino.insert(0,"?")
# print(amino)
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
    for j in range(len(i[0]) - 8):
        for k in range(1, 9):
            fsrm[i[1][j]][i[0][j]][i[0][k + j]][k] += 1




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


f = open("./txtfile/dssp_predict1.txt", 'w')

predict = ""
right = ""

for i in list2:
    # leave one out procedure

    for j in range(len(i[0]) - 8):
        for k in range(1, 9):
            fsrm[i[1][j]][i[0][j]][i[0][k + j]][k] -= 1


    f.write(i[2] + '\n')
    f.write(i[0] + '\n')
    # write the right structure
    f.write(i[1] + '\n')
    right += i[1]

    # write the predict structure
    temp = gor3(i[0])
    predict += temp
    f.write(temp + '\n')

    # restore

    for j in range(len(i[0]) - 8):
        for k in range(1, 9):
            fsrm[i[1][j]][i[0][j]][i[0][k + j]][k] += 1


f.close()

enum = 0
hnum = 0
cnum = 0
for i in right:
    if i == 'C':
        cnum += 1
    elif i == 'E':
        enum += 1
    else:
        hnum += 1

ewrong = 0
hwrong = 0
cwrong = 0
num = 0
for i in range(len(predict)):
    if predict[i] == right[i]:
        num += 1
    elif right[i] == 'E':
        ewrong += 1
    elif right[i] == 'C':
        cwrong += 1
    else:
        hwrong += 1

print(num / len(predict))


def mcc(predict, right, stru):
    TP, FP, FN, TN = 0, 0, 0, 0
    for i in range(len(predict)):
        if predict[i] == stru:
            if right[i] == stru:
                TP += 1
            else:
                FP += 1
        else:
            if right[i] == stru:
                FN += 1
            else:
                TN += 1
    MCC = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    TRP = TP / (TP + FN)
    SPC = TN / (FP + TN)

    print(TP, FP, TN, FN)

    return MCC, TRP, SPC


# predict1 = 'HHHHCCCCCHHHHCCCHHCCCCHHHH'
# right1 = 'HHHHHCCCCEEEECCCEEECCCHHHH'

#print(mcc(predict, right, 'H'))
#print(mcc(predict,right,'E'))
#print(mcc(predict,right,'C'))

t2 = time.time()

# calculate the overall

#print(t2 - t1)



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



t1 = time.time()
stride = pre_process("./txtfile/stride_info.txt", "./txtfile/stride.txt")

fsr, fs, fr = fSR(stride)
# print(fs)
# print(fr)

total = fs['C'] + fs['E'] + fs['C']

# do some pre_process, read the protein sequence and structure in to a 2d list
f = open("./txtfile/stride_protein")

list = f.readlines()
list2 = []
for i in range(0, len(list), 3):
    temp = []
    temp.append(list[i + 1].strip())
    temp.append(list[i + 2].strip())
    temp.append(list[i].strip())
    list2.append(temp)
f.close()
# print(list2)
# dealing the statistics
fsrm = {}
# initialize
amino = []
for i in amino_dict.keys():
    amino.append(amino_dict[i])
amino.insert(0,"?")
# print(amino)
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
    for j in range(len(i[0]) - 8):
        for k in range(1, 9):
            fsrm[i[1][j]][i[0][j]][i[0][k + j]][k] += 1




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


f = open("./txtfile/stride_predict1.txt", 'w')

predict = ""
right = ""

for i in list2:
    # leave one out procedure

    for j in range(len(i[0]) - 8):
        for k in range(1, 9):
            fsrm[i[1][j]][i[0][j]][i[0][k + j]][k] -= 1


    f.write(i[2] + '\n')
    f.write(i[0] + '\n')
    # write the right structure
    f.write(i[1] + '\n')
    right += i[1]

    # write the predict structure
    temp = gor3(i[0])
    predict += temp
    f.write(temp + '\n')

    # restore

    for j in range(len(i[0]) - 8):
        for k in range(1, 9):
            fsrm[i[1][j]][i[0][j]][i[0][k + j]][k] += 1


f.close()

enum = 0
hnum = 0
cnum = 0
for i in right:
    if i == 'C':
        cnum += 1
    elif i == 'E':
        enum += 1
    else:
        hnum += 1

ewrong = 0
hwrong = 0
cwrong = 0
num = 0
for i in range(len(predict)):
    if predict[i] == right[i]:
        num += 1
    elif right[i] == 'E':
        ewrong += 1
    elif right[i] == 'C':
        cwrong += 1
    else:
        hwrong += 1

print(num / len(predict))


def mcc(predict, right, stru):
    TP, FP, FN, TN = 0, 0, 0, 0
    for i in range(len(predict)):
        if predict[i] == stru:
            if right[i] == stru:
                TP += 1
            else:
                FP += 1
        else:
            if right[i] == stru:
                FN += 1
            else:
                TN += 1
    MCC = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    TRP = TP / (TP + FN)
    SPC = TN / (FP + TN)

    print(TP, FP, TN, FN)

    return MCC, TRP, SPC

# predict1 = 'HHHHCCCCCHHHHCCCHHCCCCHHHH'
# right1 = 'HHHHHCCCCEEEECCCEEECCCHHHH'

#print(mcc(predict, right, 'H'))
#print(mcc(predict,right,'E'))
#print(mcc(predict,right,'C'))

t2 = time.time()

# calculate the overall

#print(t2 - t1)





