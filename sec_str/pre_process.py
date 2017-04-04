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


f = open("amino_dict.txt")
amino_dict = {}
for i in f.readlines():
    temp = i.strip().split()
    amino_dict[temp[0].upper()] = temp[1]
f.close()

f = open("stru_dict.txt")
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
    fw.close()

    return list_2





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
    for i in list_2:
        # calculate fSR
        if i[4] in result.keys():
            if i[3] in result[i[4]].keys():
                result[i[4]][i[3]] += 1
            else:

                result[i[4]][i[3]] = 1
        else:
            temp = {}
            temp[i[3]] = 1
            result[i[4]] = temp
        # calcuate fS
        if i[4] in fS.keys():
            fS[i[4]] += 1
        else:
            fS[i[4]] = 1
        # calculate fR
        if i[3] in fR.keys():
            fR[i[3]] += 1
        else:
            fR[i[3]] = 1

    return result, fS, fR


