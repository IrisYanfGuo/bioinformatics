from numpy import sqrt

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
                if predict[i]==right[i]:
                    TN += 1
    print(TP,FP,FN,TN)
    MCC = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    TRP = TP / (TP + FN)
    SPC = TN / (FP + TN)

    return MCC

mcc('CHHHCCCCEEEECCCEEECCCHHHHC','HHHHHCCCCEEEECCCEEECCCHHHH','E')