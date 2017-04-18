from pre_process import *
from numpy import log
import matplotlib.pyplot as plt

dssp = pre_process("./txtfile/dssp_info.txt", "./txtfile/dssp_protein.txt")
stride = pre_process("./txtfile/stride_info.txt", "./txtfile/stride_protein.txt")

fsr,fs,fr= fSR(dssp)
#print(fs)
#print(fr)
total = fs['C']+fs['E']+fs['C']

# begin gor3 algorithm


def gor2(dssp,resultfile):
    predict2 = []
    for i in range(len(dssp)):

        helix = log(fsr['H'][dssp[i][3]] / (fsr['C'][dssp[i][3]] + fsr['E'][dssp[i][3]])) + log(
            (fs['E'] + fs['C']) / fs['H'])
        for j in range(-8, 9):
            if j != 0:
                if i + j > 0 and i + j < len(dssp) - 1:
                    # here now use the gor ii
                    t = i + j
                    if dssp[t][0] == dssp[i][0] and dssp[t][1] == dssp[i][1]:
                        helix += log(fsr['H'][dssp[t][3]] / (fsr['C'][dssp[t][3]] + fsr['E'][dssp[t][3]])) + log(
                            (fs['E'] + fs['C']) / fs['H'])

        coil = log(fsr['C'][dssp[i][3]] / (fsr['H'][dssp[i][3]] + fsr['E'][dssp[i][3]])) + log(
            (fs['E'] + fs['H']) / fs['C'])
        for j in range(-8, 9):
            if j != 0:
                if i + j > 0 and i + j < len(dssp) - 1:
                    # here now use the gor ii
                    t = i + j
                    if dssp[t][0] == dssp[i][0] and dssp[t][1] == dssp[i][1]:
                        coil = log(fsr['C'][dssp[t][3]] / (fsr['H'][dssp[t][3]] + fsr['E'][dssp[t][3]])) + log(
                            (fs['E'] + fs['H']) / fs['C'])

        sheet = log(fsr['E'][dssp[i][3]] / (fsr['H'][dssp[i][3]] + fsr['C'][dssp[i][3]])) + log(
            (fs['C'] + fs['H']) / fs['E'])
        for j in range(-8, 9):
            if j != 0:
                if i + j > 0 and i + j < len(dssp) - 1:
                    # here now use the gor ii
                    t = i + j
                    if dssp[t][0] == dssp[i][0] and dssp[t][1] == dssp[i][1]:
                        coil = log(fsr['E'][dssp[t][3]] / (fsr['H'][dssp[t][3]] + fsr['C'][dssp[t][3]])) + log(
                            (fs['C'] + fs['H']) / fs['E'])

        if coil == max(sheet, coil, helix):
            predict2.append('C')
        elif sheet == max(sheet, coil, helix):
            predict2.append('E')

        else:

            predict2.append('H')

    f = open(resultfile, 'w')
    f.write("".join(predict2))
    f.write("\n")

    right = ''.join([dssp[i][4] for i in range(len(dssp))])
    print(''.join(predict2))
    f.write(right)
    f.close()
    print(right)




def gor2_list(alist):
    predict2 = []
    for i in range(len(alist)):

        helix = log(fsr['H'][alist[i]] / (fsr['C'][alist[i]] + fsr['E'][alist[i]])) + log(
            (fs['E'] + fs['C']) / fs['H'])
        for j in range(-8, 9):
            if j != 0:
                if i + j > 0 and i + j < len(alist) - 1:
                    # here now use the gor ii
                    t = i + j

                    helix += log(fsr['H'][alist[t]] / (fsr['C'][alist[t]] + fsr['E'][alist[t]])) + log(
                        (fs['E'] + fs['C']) / fs['H'])

        coil = log(fsr['C'][alist[i]] / (fsr['H'][alist[i]] + fsr['E'][alist[i]])) + log(
            (fs['E'] + fs['H']) / fs['C'])
        for j in range(-8, 9):
            if j != 0:
                if i + j > 0 and i + j < len(alist) - 1:
                    # here now use the gor ii
                    t = i + j

                    coil = log(fsr['C'][alist[t]] / (fsr['H'][alist[t]] + fsr['E'][alist[t]])) + log(
                        (fs['E'] + fs['H']) / fs['C'])

        sheet = log(fsr['E'][alist[i]] / (fsr['H'][alist[i]] + fsr['C'][alist[i]])) + log(
            (fs['C'] + fs['H']) / fs['E'])
        for j in range(-8, 9):
            if j != 0:
                if i + j > 0 and i + j < len(alist) - 1:
                    # here now use the gor ii
                    t = i + j

                    coil = log(fsr['E'][alist[t]] / (fsr['H'][alist[t]] + fsr['C'][alist[t]])) + log(
                        (fs['C'] + fs['H']) / fs['E'])

        if coil == max(sheet, coil, helix):
            predict2.append('C')
        elif sheet == max(sheet, coil, helix):
            predict2.append('E')

        else:

            predict2.append('H')
    return ''.join(predict2)

#gor2(dssp,'predict2.txt')


def printresult(str,str2):
    plt.figure(figsize=(len(str),2))
    for i in range(len(str)):
        if str[i]=='C':
            plt.plot(i, 1, label='o', marker="+", markersize=30, color='r')
        elif str[i]=='E':
            plt.plot(i, 1, label='o', marker="+", markersize=30, color='g')
        else:
            plt.plot(i, 1, label='o', marker="+", markersize=30, color='b')
    for i in range(len(str2)):
        if str2[i]=='C':
            plt.plot(i, 2, label='o', marker="+", markersize=30, color='r')
        elif str2[i]=='E':
            plt.plot(i, 2, label='o', marker="+", markersize=30, color='g')
        else:
            plt.plot(i, 2, label='o', marker="+", markersize=30, color='b')
    plt.show()

#printresult(gor2_list(alistt),a)

def printwrong(str,str2):
    plt.figure(figsize=(len(str2),2))
    rate=0
    for i in range(len(str)):
        if str[i]==str2[i]:
            plt.plot(i, 2, label='o', marker="+", markersize=30, color='b')
            rate +=1
        else:

            plt.plot(i, 2, label='o', marker="+", markersize=30, color='r')
    print("right",rate/len(str) )
    plt.show()

