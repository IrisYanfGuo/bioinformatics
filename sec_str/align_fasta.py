from read_fasta import read_fasta


seq = read_fasta("./fasta/a_align.fasta")

for i in seq:
    print(i)
