inFile = open('TCGA-Cosmic')
ouFile = open('TCGA-Cosmic-Sample-Number', 'w')
D = {}
for line in inFile:
    line = line.strip()
    fields = line.split('\t')
    cancer_type = fields[7]
    sample_id = fields[4]
    D.setdefault(cancer_type, [])
    D[cancer_type].append(sample_id)
inFile.close()

N = 0
for k in D:
    ouFile.write(k +'\t' + str(len(set(D[k]))) + '\n')
    N += len(set(D[k]))
    ouFile.write('total:\t' + str(N) + '\n')
ouFile.close()
