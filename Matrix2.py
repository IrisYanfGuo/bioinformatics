# bioinformatics
# implementing the sequence alignment algorithms
# the needleman-wunsch(global) and smith-waterman(local) alighment algorithms
# dynamic programming and backtracking

# this matrix can only be used to represent the BLOSUM50 and something like that.
class Matrix2(object):
    def __init__(self,dict):
        self.matrix = dict


    def get(self,row,col):
        return self.matrix[row][col]

    def print(self):
        for key in self.matrix:
            for key2 in self.matrix[key]:

                print("{0:3d}".format(self.matrix[key][key2]),end="")

            print()



#A = Matrix("BLOSUM50.txt")
#A.print()