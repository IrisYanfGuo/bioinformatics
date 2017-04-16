from read_fasta import read_fasta
from gor3 import gor3

d = read_fasta("./fasta/d.fasta")

for i in d:
    print(i)