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



# print(len(dssp),len(stride),len(cath))


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








def pre_process(file, newfile):
    '''
    do preprocessing, only for dssp_info and stride_info
    do with uncertain values
    :param list_2: 2-d list
    :return: 1-d list
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


fsrm = {}
# initialize
amino = []
for i in amino_dict.keys():
    amino.append(amino_dict[i])
amino.insert(0,"?")


# print(stride)

def fSR(list_2):
    '''
    calculate the fSR used in the GOR III Aalgorithm
    :param list_2:  list from the reading of stride or dssp
    :return: a 2 d dictionary
    '''
    result = {}
    fS = {}
    fR = {}

    # initialize result
    for i in ['C','E','H']:
        temp ={}
        for j in amino:
            temp[j] =0.001
        result[i] = temp

    # initialize fs

    for i in ['C', 'E', 'H']:
        fS[i] =0
    for i in amino:
        fR[i] = 0

    for i in list_2:
        # calculate fSR
        result[i[4]][i[3]] +=1
        fS[i[4]] +=1
        fR[i[3]] +=1

    return result, fS, fR

def fSR2(list_2_2):
    result = {}
    fS = {}
    fR = {}
    for i in list_2_2:
        # calculate fSR
        if i[1] in result.keys():
            if i[0] in result[i[1]].keys():
                result[i[1]][i[0]] += 1
            else:

                result[i[1]][i[0]] = 1
        else:
            temp = {}
            temp[i[0]] = 1
            result[i[1]] = temp
        # calcuate fS
        if i[1] in fS.keys():
            fS[i[1]] += 1
        else:
            fS[i[1]] = 1
        # calculate fR
        if i[0] in fR.keys():
            fR[i[0]] += 1
        else:
            fR[i[0]] = 1

    return result, fS, fR



