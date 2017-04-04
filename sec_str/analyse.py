from numpy import sqrt

f = open("predict2.txt")
predict = f.readline().strip()
right = f.readline()
print(predict)
print(right)

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

print(mcc(predict, right, 'H'))
print(mcc(predict,right,'E'))
print(mcc(predict,right,'C'))
print((8 * 11 - 6 * 1) / sqrt((14) * (9) * (17) * 12))
