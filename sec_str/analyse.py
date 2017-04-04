f = open("predict1.txt")
predict = f.readline()
right = f.readline()
print(predict)
print(len(predict))
print(len(right))
print(right)

enum = 0
hnum =0
cnum = 0
for i in right:
    if i == 'C':
        cnum +=1
    elif i == 'E':
        enum +=1
    else:
        hnum +=1


ewrong = 0
hwrong = 0
cwrong = 0
num = 0
for i in range(len(predict)-1):
    if predict[i] == right[i]:
        num += 1
    elif right[i] == 'E':
        ewrong +=1
    elif right[i] == 'C':
        cwrong +=1
    else:
        hwrong +=1

print(num/len(predict))
print(ewrong/enum)
print(cwrong/cnum)
print(hwrong/hnum)
