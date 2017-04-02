def read_fasta(filename):
    f = open(filename)
    proteins = []
    d = f.read().split(">")

    # deal with the  "" before the first ""
    d.pop(0)
    for i in d:
        t = i.splitlines()
        # get rid of the line of the protein name
        t.pop(0)
        s = "".join(t)
        proteins.append(s)
    return proteins

# test
# j = read_fasta("WW-homo-136.fasta")
# print(len(j),j[0])
# print(len(j[0]))
# print(read_fasta("protein-sequences.fasta"))
