inFile = open('TCGA-Cosmic-Count')
ouFile = open('TCGA-Cosmic-Count-mutation-type', 'w')
D = {}
for line in inFile:
    line = line.strip()
    fields = line.split('\t')
    D[fields[3]] = 1
inFile.close()
for x in sorted(D):
    ouFile.write(x + '\n')
ouFile.close()


inFile = open('TCGA-Cosmic-Count')
for line in inFile:
    line = line.strip()
    fields = line.split('\t')
    if fields[3] == '':
        print(line)
inFile.close()


