L = []
ouFile = open('GWAS-Gene', 'w')
inFile = open('gwascatalog-2015.1.28.txt')
head = inFile.readline()
for line in inFile:
    line = line.strip()
    fields = line.split('\t')
    gene = fields[14].strip()
    if gene == '-':
        pass
    elif gene.find(';') != -1:
        gs = gene.split(';')
        for g in gs:
            if g.find(' - ') != -1:
                t = g.split(' - ')
                L+=t
            else:
                L.append(g)
inFile.close()
S = set(L)
for x in S:
    ouFile.write(x + '\n')
ouFile.close()
