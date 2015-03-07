inFile = open('TCGA-Cosmic')
ouFile = open('TCGA-Sample-NonStopGain', 'w')
ouFile2 = open('TCGA-Sample-NonStopGain-Cancer', 'w')
D = {}
for line in inFile:
    line = line.strip()
    fields = line.split('\t')
    gene = fields[0].split('_')[0]
    sample = fields[4]
    cancer = fields[7]
    mutation = fields[15]
    #if mutation == 'Substitution - Nonsense' or mutation == 'Nonstop extension':
    if mutation == 'Substitution - Nonsense':
        D.setdefault(cancer, {})
        D[cancer].setdefault(sample, [])
        D[cancer][sample].append(gene)
inFile.close()

for cancer in D:
    d = D[cancer].items()
    ouFile2.write(cancer + '\t' + str(len(D[cancer])) + '\n')
    d.sort(cmp = lambda x,y:cmp(len(set(x[1])), len(set(y[1]))), reverse = True)
    for item in d:
        ouFile.write(cancer + '\t' + item[0] + '\t' + str(len(set(item[1]))) +'\t' + '\t'.join(set(item[1])) + '\n')



#d = D.items()
#d.sort(cmp = lambda x,y:cmp(len(set(x[1])), len(set(y[1]))), reverse = True)

#for item in d:
#    ouFile.write(item[0] + '\t' + str(len(set(item[1]))) +'\t' + '\t'.join(set(item[1])) + '\n')

ouFile.close()
ouFile2.close()

