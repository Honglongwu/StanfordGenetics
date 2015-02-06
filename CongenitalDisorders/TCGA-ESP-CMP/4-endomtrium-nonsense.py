inFile = open('TCGA-ESP-count-enriched-significant')
ouFile = open('TCGA-ESP-count-enriched-significant-endometrium-nonsense', 'w')
for line in inFile:
    line = line.strip()
    fields = line.split('\t')
    mutation = fields[6]
    cancer = fields[2]
    if mutation == 'Nonsense' and cancer == 'endometrium':
        ouFile.write(line + '\n')

inFile.close()
ouFile.close()
