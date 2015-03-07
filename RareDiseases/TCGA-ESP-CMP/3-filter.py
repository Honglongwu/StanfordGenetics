inFile = open('TCGA-ESP-count-enriched-total')
ouFile = open('TCGA-ESP-count-enriched-total.filtered', 'w')
for line in inFile:
    line = line.strip()
    fields = line.split('\t')
    gene = fields[0]
    if gene.find('_')==-1:
        ouFile.write(line + '\n')
inFile.close()
ouFile.close()
