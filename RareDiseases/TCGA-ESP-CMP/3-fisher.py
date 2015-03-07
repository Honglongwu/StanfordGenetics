import scipy.stats
Sig = 0.05
inFile = open('TCGA-ESP-count-enriched')
ouFile = open('TCGA-ESP-count-enriched-significant', 'w')
ouFile2 = open('TCGA-ESP-count-enriched-total', 'w')
L = []
L2 = []
for line in inFile:
    line = line.strip()
    fields = line.split('\t')
    try:
        n1 = int(fields[4])
        n2 = int(fields[5])
        n3 = int(fields[10])
        n4 = int(fields[11])
        p = scipy.stats.fisher_exact([[n1,n2],[n3,n4]])[1]
        if p < Sig:
            L.append(fields + [p])
        L2.append(fields + [p])
    except:
        pass
inFile.close()

L.sort(cmp = lambda x,y:cmp(x[-1], y[-1]))
L2.sort(cmp = lambda x,y:cmp(x[-1], y[-1]))
for x in L:
    ouFile.write('\t'.join([str(y) for y in x]) + '\n')
for x in L2:
    ouFile2.write('\t'.join([str(y) for y in x]) + '\n')

ouFile.close()
ouFile2.close()
