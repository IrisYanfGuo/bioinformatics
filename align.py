# real things
from Matrix import *
from read_fasta import *

seq_list = read_fasta("WW-sequence.fasta")
sequ1 = seq_list[0]
sequ2 = seq_list[1]
mat1 = Matrix("BLOSUM50.txt")
print(sequ1, sequ2)


# mat1.print()





def score(seq1, seq2, mat):
    if (len(seq1) == 0 and len(seq2) == 0):
        return 0
    elif len(seq1) == 0:
        return score(seq1, seq2[:-1], mat) + mat.get(seq2[-1], '*')
    elif len(seq2) == 0:
        return score(seq1[:-1], seq2, mat) + mat.get(seq1[-1], '*')
    else:
        return max(score(seq1, seq2[:-1], mat) + mat.get(seq2[-1], '*'),
                   score(seq1[:-1], seq2, mat) + mat.get(seq1[-1], '*'),
                   score(seq1[:-1], seq2[:-1], mat) + mat.get(seq1[-1], seq2[-1]))


print(score(sequ1, sequ2, mat1))
