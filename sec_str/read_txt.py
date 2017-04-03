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


cath = read_txt("cath_info.txt")
stride = read_txt("stride_info.txt")
dssp = read_txt("dssp_info.txt")

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


def pre_process(list_2):
    '''
    do preprocessing, change the 3 character amino to 1 character
    do with uncertain values
    :param list_2: 2-d list
    :return: 1-d list
    '''
    pass
