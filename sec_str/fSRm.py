
# do some pre_process, read the protein sequence and structure in to a 2d list
f = open("dssp_protein")

list = f.readlines()
list2 =[]
for i in range(0,len(list),3):
    temp =[]
    temp.append(list[i+1].strip())
    temp.append(list[i+2].strip())
    list2.append(temp)

# dealing the statistics
fsrm={}
for i in list2:
    for j in range(len(i[0])):
        for k in range(1,9):
            if i[1][j] in fsrm.keys():
                if i[0][j] in fsrm[i[1][j]].keys():
                    if i[0][k] in fsrm[i[1][j]][i[0][j]].keys():
                        fsrm[i[1][j]][i[0][j]][i[0][k]] +=1
                    else:

                        temp =1
                        fsrm[i[1][j]][i[0][j]][i[0][k]] = temp
                else:

                    temp = 1
                    temp2 ={}
                    temp2[i[0][k]] =temp
                    fsrm[i[1][j]][i[0][j]] = temp2
            else:

                temp= 1
                temp2 = {}
                temp2[i[0][j]] = temp
                temp3 ={}
                temp3[i[0][j]] = temp2

                fsrm[i[1][j]] = temp3

for i in fsrm.keys():
    for j in fsrm[i].keys():
        print(fsrm[i][j])









