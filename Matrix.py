# bioinformatics
# implementing the sequence alignment algorithms
# the needleman-wunsch(global) and smith-waterman(local) alighment algorithms
# dynamic programming and backtracking


class Matrix(object):
    def __init__(self,file):
        f = open(file)
        self.matrix ={}
        # c is the name of the column
        c = f.readline().split()
        for i in f:
            r = i.split()
            rowname = r.pop(0)
            t = {}
            for j in range(len(c)):
                t[c[j]] = r[j]
            self.matrix[rowname] = t


    def __str__(self):
        s =''
        for key in self.matrix:
            pass


    def get(self,row,col):
        return self.data[row][col]



A = Matrix("BLOSUM50.txt")
print(A)
