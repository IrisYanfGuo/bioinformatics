# bioinformatics
# implementing the sequence alignment algorithms
# the needleman-wunsch(global) and smith-waterman(local) alighment algorithms
# dynamic programming and backtracking


class Matrix(object):
    def __init__(self,row,col):
        self.data =[[0]*col for i in range(row)]

    def __str__(self):
        s =''
        for i in range(len(self.data)):
            s += str(self.data[i])+'\n'
        return s

    def set(self,row,col,n):
        self.data[row][col] =n

    def get(self,row,col):
        return self.data[row][col]





f = open("BLOSUM50.txt")
a = Matrix(2,3)
for i in f:
    if(i[0]!='#'):
        print(i.strip().split(" "))
        print(len(i.strip().split("  ")))




