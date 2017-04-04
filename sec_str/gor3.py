from pre_process import *
from numpy import log

dssp = pre_process("dssp_info.txt", "dssp.txt")
stride = pre_process("stride_info.txt", "stride.txt")

fsr,fs,fr= fSR(dssp)
print(fs)
print(fr)
total = fs['C']+fs['E']+fs['C']

# begin gor3 algorithm


def gor2(dssp,resultfile):
    predict2 = []
    for i in range(len(dssp)):

        helix = log(fsr['H'][dssp[i][3]] / (fsr['C'][dssp[i][3]] + fsr['E'][dssp[i][3]])) + log(
            (fs['E'] + fs['C']) / fs['H'])
        for j in range(-5, 6):
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

gor2(stride,'predict2.txt')


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


alistt = 'ITKVEAENMKIGGTYAGKISAPFDGVALYANADYVSYSQYFANSTHNISVRGASSNAGTAKVDLVIGGVTVGSFNFTGKTPTVQTLSNITHATGDQEIKLALTSDDGTWDAYVDFIEFSL'
print(gor2_list(alistt))
print('CEEEECCCCEEECCCCEEECCCCCEEEECCCCCEEEEEEEECCCEEEEEEEEEECCCCEEEEEEEECCEEEEEEEEECCCCEEEEEEEEECCCEEEEEEEEECCCCCCCCEEEEEEEEEC')