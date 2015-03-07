D = {}
inFile = open('ESP-snp-indels-Count-unified')
for line in inFile:
    line = line.strip()
    fields = line.split('\t')
    gene = fields[0]
    D.setdefault(gene, [])
    D[gene].append(fields)
inFile.close()

inFile = open('TCGA-Cosmic-Count-unified')
ouFile = open('TCGA-ESP-count', 'w')
for line in inFile:
    line = line.strip()
    fields = line.split('\t')
    mutation_type = fields[-1]
    gene = fields[0].split('_')[0]
    enst = fields[1]
    if gene in D:
        for x in D[gene]:
            if mutation_type == x[-2] and enst == x[-1]:
                ouFile.write(line + '\t' + '\t'.join(x) + '\n')

inFile.close()
ouFile.close()
