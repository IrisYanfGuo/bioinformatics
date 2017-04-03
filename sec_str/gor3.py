from pre_process import *
from numpy import log

dssp = pre_process("dssp_info.txt", "dssp.txt")
stride = pre_process("stride_info.txt", "stride.txt")

fsr,fs,fr= fSR(dssp)
print(fs)
print(fr)
total = fs['C']+fs['E']+fs['C']

# begin gor3 algorithm

predict2 =[]
for i in range(len(dssp)):

    helix = log(fsr['H'][dssp[i][3]]/(fsr['C'][dssp[i][3]]+fsr['E'][dssp[i][3]])) + log((fs['E']+fs['C'])/fs['H'])
    for j in range (-8,9):
        if j != 0:
            if i+j > 0 and i+j < len(dssp)-1:
                # here now use the gor ii
                t = i+j
                if dssp[t][0] == dssp[i][0] and dssp[t][1] ==dssp[i][1]:
                    helix += log(fsr['H'][dssp[t][3]]/(fsr['C'][dssp[t][3]]+fsr['E'][dssp[t][3]]))+log((fs['E']+fs['C'])/fs['H'])

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


    if coil == max(sheet,coil,helix):
        predict2.append('C')
    elif sheet == max(sheet,coil,helix):
        predict2.append('E')

    else:

        predict2.append('H')

f = open("predict1.txt",'w')
f.write("".join(predict2))
f.write("\n")



right = ''.join([dssp[i][4] for i in range(len(dssp))])
print(''.join(predict2))
f.write(right)
f.close()
print(right)



