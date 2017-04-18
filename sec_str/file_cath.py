import matplotlib.pyplot as  plt
from numpy import sqrt

f = open("./txtfile/cath_info.txt")

fw = open("./txtfile/cath.txt", 'w')

dict = {}
dict['Alpha'] = 'A'
dict['Beta'] = 'B'
dict["Alpha/beta"] = 'AB'
dict["None"] = 'N'

for i in f.readlines():
    t = i.strip().split()

    t[2] = dict[t[2]]

    fw.write(t[0] + ' ')
    fw.write(t[1] + ' ')
    fw.write(t[2] + '\n')

fw.close()
f.close()

f = open("./txtfile/cath.txt")
right = []
for i in f.readlines():
    t = []
    t.append(i.strip().split()[0])
    t.append(i.strip().split()[2])
    right.append(t)

f.close()

f = open("./txtfile/dssp_predict.txt")
list = f.readlines()
list2 = []
for i in range(0, len(list), 4):
    t = []
    t.append(list[i].strip().split()[0])
    t.append(list[i].strip().split()[2])
    t.append(list[i + 3].strip())
    t.append(list[i + 2].strip())
    t.append(list[i + 1].strip())
    list2.append(t)


# print(list2)

def family_predict(alist):
    cnum = 1
    enum = 1
    hnum = 1
    for i in alist:
        if i == 'C':
            cnum += 1
        elif i == 'E':
            enum += 1
        else:
            hnum += 1

    return (hnum-enum)/(len(alist)-cnum)






result = []
for i in list2:
    t = []
    t.append(family_predict(i[2]))
    t.append(i[0])
    t.append(i[2])
    t.append(i[3])
    t.append(i[4])
    result.append(t)

#print(result)


def q3(right, predict):
    num = 0
    for i in range(len(right)):
        if right[i] == predict[i]:
            num += 1
    return num / len(right)


def mcc(predict, right, stru):
    TP, FP, FN, TN = 1, 1, 1, 1
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


temp = []
accuracy = 0
for i in range(len(right)):
    t = []
    t.append(result[i][0])
    t.append(right[i][1])

    if result[i][0] > 0.45:
        t.append('A')
    elif result[i][0] < -0.25:
        t.append('B')
    else:
        t.append('AB')

    t.append(right[i][0])

    t.append(result[i][2])
    t.append(result[i][3])
    t.append(result[i][4])
    t.append(q3(result[i][3], result[i][2]))
    t.append(mcc(result[i][3], result[i][2], 'H'))
    t.append(mcc(result[i][3], result[i][2], 'E'))
    t.append(mcc(result[i][3], result[i][2], 'C'))
    temp.append(t)

    if t[1] == t[2]:
        accuracy += 1

print("accuracy = ", accuracy / len(right))

for i in temp:
    if i[1] == 'A':
        plt.plot([1], [i[0]], 'ro')
    elif i[1] == 'B':
        plt.plot([2], [i[0]], 'bo')
    else:
        plt.plot([3], [i[0]], 'go')

# plt.show()





temp.sort()
f = open("./txtfile/cathinfo2.csv", 'w')
f.write('H-E,')
f.write('actual family,')
f.write('predict family,')
f.write('name,')
f.write('predict structure,')
f.write('actual structure,')
f.write("amino,")
f.write('q3,')
f.write('mccH,')
f.write('mccE,')
f.write('mccC,')
f.write('\n')

for i in temp:
    f.write(str(i[0]) + ',')
    f.write(str(i[1]) + ',')
    f.write(str(i[2]) + ',')
    f.write(str(i[3]) + ',')
    f.write(str(i[4]) + ',')
    f.write(str(i[5]) + ',')
    f.write(str(i[6]) + ',')
    f.write(str(i[7]) + ',')
    f.write(str(i[8]) + ',')
    f.write(str(i[9]) + ',')
    f.write(str(i[10]) + '\n')

#print(temp)
