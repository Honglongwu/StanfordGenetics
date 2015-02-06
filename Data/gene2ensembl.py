inFile = open('gene2ensembl')
ouFile = open('gene2ensembl-human', 'w')
for line in inFile:
    line = line.strip()
    fields = line.split('\t')
    if fields[0] == '9606':
        ouFile.write(line + '\n')
inFile.close()
ouFile.close()
