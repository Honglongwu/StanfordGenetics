#inFile = open('TCGA-Cosmic')
inFile = open('TCGA-Cosmic-part')
ouFile = open('TCGA-Cosmic-endometrium-nonsense', 'w')
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

'''
for cancer in D:
    for sample in D[cancer]:
        print('\t'.join(set(D[cancer][sample])))
'''

#d = D.items()
#d.sort(cmp = lambda x,y:cmp(len(set(x[1])), len(set(y[1]))), reverse = True)

#for item in d:
#    ouFile.write(item[0] + '\t' + str(len(set(item[1]))) +'\t' + '\t'.join(set(item[1])) + '\n')

ouFile.close()

