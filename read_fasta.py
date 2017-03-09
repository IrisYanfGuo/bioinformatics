f = open("WW-sequence.fasta")
d =[]
for i in f:
    i

#print(d)
f.close()

f2= open("WW-sequence.fasta")
proteins = []
d = f2.read().split(">")
d.pop(0)
for i in d:
    t = i.splitlines()
    t.pop(0)
    s = ''.join(t)
    proteins.append(s)

for i in range(len(proteins)):
    print(proteins[i])