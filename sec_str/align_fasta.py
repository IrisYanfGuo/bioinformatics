from read_fasta import read_fasta
from gor3 import gor3




# GHR--G-D -> GHRGD
def tostr(alist):
    result = ''
    for i in alist:
        if i != '-':
            result += i
    return result


# GHRGD -->EEE--C-D
def reverse(alist1, alist):
    result = ''

    # index of alist
    j = 0
    for i in alist1:
        if i != '-':
            result += alist[j]
            j += 1
        else:
            result += '-'
    return result


# print(reverse(seq[1][start:end],tostr(seq[1][start:end])))

# print(seq[1][start:end])
# print(tostr(seq[1][start:end]))





astru_right="CCCCCCCCCCCCCCHHHHHHHHHHHHHHCCCCEEEEEEEECCCCCEEEEEEECCCCCCCCEEEEEECCCCCCHHHHHHHHHHHHHHHHHCCCCHHHHHHHHHCEEEEECCCCHHHHHHHHHCCCCCCCCCCCCCCCCCCCCCHHHCCCCCCCCCCCECCCCCCCECCCCCCCCHHHHHHHHHHHHHCCEEEEEEEEECCCEEEECCCCCCCCCCCHHHHHHHHHHHHHHHHHHHCCCCEEEEHHHHCCCCCCCHHHHHHHCCCCEEEEEEECCCCCCHHHCCHHHHHHHHHHHHHHHHHHHHHHHHC"
bstru_right ='CCCCECECCCCCECECCCEEEEEECCHHHCCCEEEEEECCEEEEEEEEECCCCCCCCCCEEEEECCCCCCCCECECCCCEEEEECCCCCCCCCCECEECCCCECCECEEECCCCCCCCCCCHHHCEEEEECECCCCCCEEEEEECCCEEECEEECCCCCCCCEEECCCCECCEEEEEECCC'
cstru_right = 'CCCCCCCCCCCCCCCCCCCCHHHHHHHHHCHHHHHHHHHHHHHHHHHHHHHHCCCCCHHHHHHHHHHHHHHHHHHHHHHHCECCCCCCCCCCCHHHHHHHHHHHCCHHHHHHHHHHHHHCCCCCEEEEEEEECCCCEEEEEEEECCCECCCCCCEEEEEEECCHHHCHHHHCCCHHHHHHHHHHHECHHHHHHHHHHHCCC'
dstru_right = 'CCCCCECCCECCCECCCCCCCEEEEEEECCEEEEEEEHHHHHHHHHHHHHHHHHHHCCCCEECCCCCCCCCCCCCHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHEEECCCCEEEECCCCCHHHHHHHHCCCCCCHHHHHHHHHHHCCCCC'
estru_right ='CCCCCCCCCCCHHHCCCCHHHHHHHHHHHHHHHHHHHCCCHHHHHHHHHHHHHHHHHHHHHHHHCCCCCHHHHHHHHHHHHHHHHCCCCCHHHHHHHHHHHHHHCHHHHCCCC'
fstru_right = 'CCCCCCCCCCCCCCCHHHHHHHHHHHHCCCCHHHHHHHCEECCCCEECCCCCCCEEEECCCCCHHHHCCCCCCCCCCCCCCCCC'
gstru_right = 'CCHHHHHHHCCCCCCCCCCCCCCCCCHHHHHHHCCCCHHHHHHHHCCCCCCCCCCCCHHHHCCCCCCCCCCCHHHHHHHHHHHHCC'
'''
f = open("./fasta/6sequence.txt",'w')
f.write(tostr(read_fasta("./fasta/a_align.fasta")[0]))
f.write("\n")
f.write(tostr(read_fasta("./fasta/b_align.fasta")[0]))
f.write("\n")
f.write(tostr(read_fasta("./fasta/c_align.fasta")[0]))
f.write("\n")
f.write(tostr(read_fasta("./fasta/d_align.fasta")[0]))
f.write("\n")

f.write(tostr(read_fasta("./fasta/e_align.fasta")[0]))
f.write("\n")

f.write(tostr(read_fasta("./fasta/f_align.fasta")[0]))
f.write("\n")
f.write(astru_right+"\n")
f.write(bstru_right+"\n")
f.write(cstru_right+"\n")
f.write(dstru_right+"\n")
f.write(estru_right+"\n")
f.write(fstru_right+"\n")
'''

def q3(list1,list2):
    right = 0
    for i in range(len(list1)):
        if list1[i] == list2[i]:
            right+=1
    return right/len(list1)




def align_predict(filename):
    seq = read_fasta(filename)
    for i in range(len(seq[0])):
        if seq[0][i] != '-':
            start = i
            break

    for i in range(len(seq[0]) - 1, 0, -1):
        if seq[0][i] != '-':
            end = i + 1
            break

    result = []
    for i in range(len(seq)):
        t1 = seq[i][start:end]
        t2 = gor3(tostr(t1))
        result.append(reverse(t1, t2))
    '''
    for i in result:
        print(i)

    '''
    stru_pred = ''
    for i in range(len(result[0])):
        if result[0][i] != '-':
            temp = [result[j][i] for j in range(len(result))]

            t1 = temp.count('E')
            t2 = temp.count('C')
            t3 = temp.count('H')

            if t1 == max(t1, t2, t3):
                stru_pred += 'E'
            elif t2 == max(t1, t2, t3):
                stru_pred += 'C'
            else:
                stru_pred += 'H'
    return stru_pred


'''
t = align_predict("./fasta/a_align1.fasta")
print(t)
print(astru_right)
print(q3(t,astru_right))
print(q3(gor3(tostr(read_fasta("./fasta/a_align1.fasta")[0])),astru_right))

'''





t = align_predict("./fasta/b_align3.fasta")
print(t)
print(bstru_right)
print(q3(t,bstru_right))
b = read_fasta("./fasta/b_align3.fasta")
print(q3(gor3(tostr(b[0])),bstru_right))
print(len(b))



'''
t = align_predict("./fasta/c_align2.fasta")
print(t)
print(cstru_right)
print(q3(t,cstru_right))
print(q3(gor3(tostr(read_fasta("./fasta/c_align1.fasta")[0])),cstru_right))

'''


'''
t = align_predict("./fasta/d_align.fasta")
print(t)
print(dstru_right)
print(q3(t,dstru_right))
print(q3(gor3(tostr(read_fasta("./fasta/d_align.fasta")[0])),dstru_right))
'''

'''
t = align_predict("./fasta/e_align1.fasta")
print(t)
print(estru_right)
print(q3(t,estru_right))
b = read_fasta("./fasta/e_align1.fasta")
print(q3(gor3(tostr(b[0])),estru_right))
print(len(b))

'''


'''
t = align_predict("./fasta/f_align.fasta")
print(t)
print(fstru_right)
print(q3(t,fstru_right))
print(q3(gor3(tostr(read_fasta("./fasta/f_align.fasta")[0])),fstru_right))

'''

