def sampleMissingNumber(inF):
    Sample = []
    ouFile = open(inF.split('.tab')[0] + '-missingNumber','w')
    D = {}
    inFile = open(inF)
    head = inFile.readline()
    for line in inFile:
        line = line.strip('\n')
        fields = line.split('\t')
        sample = fields[0]
        missingN = fields.count('N/A')
        D[sample] = missingN
    d = D.items()
    d.sort(cmp = lambda x,y:cmp(x[1],y[1]))
    for item in d:
        ouFile.write(item[0] + '\t' + str(item[1]) + '\n')
        if item[1] < 5000:
            Sample.append(item[0])
    ouFile.close()
    print("Sample Number:" + '\t' + str(len(Sample)))
    return Sample
def missing(inF):
    ouFile = open(inF.split('.tab')[0] + '-noMissing.tab', 'w')
    Sample = sampleMissingNumber(inF)
    D = {}
    inFile = open(inF)
    head = inFile.readline().strip('\n').split('\t')
    for x in head[1:]:
        D.setdefault(x,[])
    for line in inFile:
        line = line.strip('\n')
        fields = line.split('\t')
        sample = fields[0]
        if sample in Sample:
            for i in range(1,len(fields)):
                D[head[i]].append(fields[i])
    inFile.close()
    D2 = {}
    for k in D:
        if 'N/A' not in D[k]:
            D2[k] = D[k]
    H = []
    for g in head:  
        if g in D2:
            H.append(g)
    ouFile.write(''+'\t'+'\t'.join(H) + '\n')
    for i in range(len(Sample)):
        L = []
        for g in head:  
            if g in D2:
                L.append(D2[g][i])
                H.append(g)
        ouFile.write(Sample[i] + '\t' + '\t'.join(L) + '\n')
    ouFile.close()
    
    print("Genotype Number:" + '\t' + str(len(D2)))
missing('Yeast-Genotype.tab')
#sampleMissingNumber('Yeast-Genotype.tab')

