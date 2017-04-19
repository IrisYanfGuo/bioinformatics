from pre_process import *
from numpy import log
import time
from numpy import sqrt
import matplotlib.pyplot as plt

t1 = time.time()
dssp2 = pre_process("./txtfile/dssp_info.txt", "./txtfile/dssp.txt")
stride = pre_process("./txtfile/stride_info.txt", "./txtfile/stride.txt")

fsr, fs, fr = fSR(dssp2)

print(fsr)
print(fs)
print(fr)
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
amino.insert(0, "?")
# print(amino)
for i in ['C', 'E', 'H']:
    temp3 = {}
    for j in amino:
        temp2 = {}
        for k in amino:
            temp = {}
            for n in range(1, 9):
                # solve the divided by 0 problem
                temp[n] = 0.001

            temp2[k] = temp

        temp3[j] = temp2
    fsrm[i] = temp3

for i in list2:
    for j in range(len(i[0]) - 8):
        for k in range(1, 9):
            fsrm[i[1][j]][i[0][j]][i[0][k + j]][k] += 1

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


f = open("./txtfile/dssp_predict.txt", 'w')
def q3(right, predict):
    num = 0
    for i in range(len(right)):
        if right[i] == predict[i]:
            num += 1
    return num / len(right)

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

    return MCC


t2 = time.time()

# calculate the overall

# print(t2 - t1)

print("MCC SCORE FOR 'H' ", mcc(predict, right, 'H'))
print("MCC SCORE FOR 'E' ", mcc(predict, right, 'E'))
print("MCC SCORE FOR 'C' ", mcc(predict, right, 'C'))

f = open("./txtfile/cath.txt")
f = open("./txtfile/cath.txt")
right1 = []
for i in f.readlines():
    t = []
    t.append(i.strip().split()[0])
    t.append(i.strip().split()[2])
    right1.append(t)
print(right1)

f.close()

f = open("./txtfile/dssp_predict.txt")
list = f.readlines()
list2 = []
for i in range(0, len(list), 4):
    t = []
    t.append(list[i].strip().split()[0])
    t.append(list[i + 3].strip())
    t.append(list[i + 2].strip())
    list2.append(t)

familyA = ''
familyA_pred = ''
familyAB_pred = ''
familyB_pred = ''
familyB = ''
familyAB = ''

for i in range(len(right1)):
    if right1[i][1] == 'A':
        familyA += list2[i][1]
        familyA_pred += list2[i][2]
    elif right1[i][1] =='B':
        familyB += list2[i][1]
        familyB_pred += list2[i][2]
    else:
        familyAB += list2[i][1]
        familyAB_pred += list2[i][2]

print(q3(familyA,familyA_pred))
print(mcc(familyA,familyA_pred,'H'))
print(mcc(familyA,familyA_pred,'E'))
print(mcc(familyA,familyA_pred,'C'))
print(q3(familyB,familyB_pred))

print(q3(familyAB,familyAB_pred))
