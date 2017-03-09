# bioinformatics
# implementing the sequence alignment algorithms
# the needleman-wunsch(global) and smith-waterman(local) alighment algorithms
# dynamic programming and backtracking

# this matrix can only be used to represent the BLOSUM50 and something like that.
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
                t[c[j]] = int(r[j])
            self.matrix[rowname] = t


    def get(self,row,col):
        return self.data[row][col]

    def print(self):
        for key in self.matrix:
            for key2 in self.matrix[key]:
                
                print("{0:3d}".format(self.matrix[key][key2]),end="")

            print()



A = Matrix("BLOSUM50.txt")
A.print()