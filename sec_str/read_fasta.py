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


