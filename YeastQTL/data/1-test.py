inFile = open('Yeast-Genotype.tab')
L = []
for line in inFile:
    line = line.rstrip('\n')
    fields = line.split('\t')
    L.append(fields)
inFile.close()

for i in range(1,len(L[0])):
    D = {}
    for x in L[1:]:
        D.setdefault(x[i], 0)
        D[x[i]] += 1
    print(D)


