inFile = open('human-ngly1-interaction-gene.txt')
D = {}
for line in inFile:
    line = line.strip()
    D[line] = 1
inFile.close()

inFile = open('UCEC-3NGLY1-36Control-Expressed-Genes')
ouFile = open('UCEC-3NGLY1-36Control-Expressed-Genes-NGLY1', 'w')
for line in inFile:
    fields = line.split()
    if fields[0] in D:
        ouFile.write(line)
inFile.close()
ouFile.close()
