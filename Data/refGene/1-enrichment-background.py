D = {}
inFile = open('refGene.2015.1.23.txt')
ouFile = open('refGene-symbol.txt', 'w')
for line in inFile:
    line = line.strip()
    fields = line.split('\t')
    transcript = fields[1]
    gene = fields[12]
    if transcript.find('NM') != -1:
        D[gene] = 1
inFile.close()
for k in D:
    ouFile.write(k + '\n')
ouFile.close()
