



# transfer dssp_predict.txt file into format like fasta

def tofasta(file,newfile):
    f = open(file)
    f2 = open(newfile,'w')
    t = f.readline().strip().split()
    f2.write(t[0] + ' ' + t[1] + ' ' + t[2])
    tlist = []
    t2list = []
    tlist.append(t[3])
    t2list.append(t[4])
    for i in f.readlines():
        temp = i.strip().split()
        if temp[0] == t[0] and temp[1] == t[1]:
            tlist.append(temp[3])
            t2list.append(temp[4])
            pos = temp[2]
        else:
            f2.write('-' + pos)
            f2.write('\n')
            f2.write(''.join(tlist))
            f2.write('\n')
            f2.write(''.join(t2list))
            f2.write('\n')
            tlist = []
            t2list = []
            t = i.strip().split()
            f2.write(t[0] + ' ' + t[1] + ' ' + t[2])
            tlist.append(t[3])
            t2list.append(t[4])

    f2.write('\n')
    f2.write(''.join(tlist))
    f2.write('\n')
    f2.write(''.join(t2list))
    f2.write('\n')

    f.close()
    f2.close()

tofasta('./txtfile/dssp.txt','./txtfile/dssp_protein')
tofasta('./txtfile/stride.txt','./txtfile/stride_protein')