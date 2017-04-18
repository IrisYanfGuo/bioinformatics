 # dealing with the input and output

 # amino dict and stru dict, use to translate file

f = open("./txtfile/amino_dict.txt")
amino_dict = {}
for i in f.readlines():
    temp = i.strip().split()
    amino_dict[temp[0].upper()] = temp[1]
f.close()

f = open("./txtfile/stru_dict.txt")
stru_dict = {}
for i in f.readlines():
    temp = i.strip().split()
    stru_dict[temp[0].upper()] = temp[1]
f.close()


# read_txt
def read_txt(filename):
    '''
    :use to read structured txt file seperated with ' '
    :param filename: the position of the file
    :return: a two demension list
    '''
    f = open(filename)
    result = []
    for i in f.readlines():
        temp = i.strip().split()
        result.append(temp)
    f.close()
    return result

# dealing with the translate

# 1st translate method dssp_info --> dssp.txt

def pre_process(file, newfile):
    '''
    do preprocessing, only for dssp_info and stride_info
    do with uncertain values
    :param list_2: 2-d list
    :return: 1*5 list in this program
    '''
    list_2 = read_txt(file)

    fw = open(newfile, 'w')
    for i in list_2:
        if i[3].upper() in amino_dict.keys():
            i[3] = amino_dict[i[3].upper()]
        else:
            i[3] = "?"
        if i[4].upper() in stru_dict.keys():
            i[4] = stru_dict[i[4].upper()]
        s = ''
        for j in i:
            s += j + ' '
        s += '\n'
        fw.write(s)

    return list_2

dssp = pre_process("./txtfile/dssp_info.txt", "./txtfile/dssp_protein.txt")
stride = pre_process("./txtfile/stride_info.txt", "./txtfile/stride_protein.txt")

# 2nd method dssp.txt --> dssp_protein.txt
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


# readfasta
def read_fasta(filename):
    f = open(filename)
    proteins = []
    d = f.read().split(">")

    # deal with the  "" before the first ""
    d.pop(0)
    for i in d:
        t = i.splitlines()
        # get rid of the line of the protein name
        t.pop(0)
        s = "".join(t)
        proteins.append(s)
    return proteins



