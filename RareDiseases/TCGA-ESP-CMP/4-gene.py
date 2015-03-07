D = {}
inFile = open('TCGA-ESP-count-enriched-significant-endometrium-nonsense')
ouFile = open('TCGA-ESP-count-enriched-significant-endometrium-nonsense-gene', 'w')
for line in inFile:
    fields = line.split('\t')
    gene = fields[0].split('_')[0]
    D[gene] = 1
inFile.close()
for k in D:
    ouFile.write(k + '\n')
ouFile.close()
