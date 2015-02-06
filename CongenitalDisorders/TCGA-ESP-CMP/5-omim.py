inFile = open('../OMIM/omim-congenital')
D = {}
for line in inFile:
    line = line.strip()
    fields = line.split('|')
    gene = fields[1].split(',')[0].strip()
    D[gene] = 1
inFile.close()

inFile = open('TCGA-ESP-count-enriched-significant')
for line in inFile:
    line = line.strip()
    fields = line.split('\t')
    gene = fields[0]
    if gene in D:
        print(gene)
inFile.close()
