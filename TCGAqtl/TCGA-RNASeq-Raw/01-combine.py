import os
Fs = os.listdir('.')
D = {}
Gene = []
Sample = []
for F in Fs:
    if F[-4:] == '.txt':
        cancer = F.split('__')[0]
        platform = F.split('__')[2].split('_')[0][8:]

        inFile = open(F)
        head1 = inFile.readline().strip('\n').split('\t')
        head2 = inFile.readline().strip('\n').split('\t')
        for line in inFile:
            line = line.strip()
        #for x in range(10):
        #    line = inFile.readline().strip()
            fields = line.split('\t')
            gene = fields[0]
            if gene not in Gene:
                Gene.append(gene)
            D.setdefault(gene, {})
            for i in range(2,len(fields),2):
                sampleID = head1[i] + '-' +cancer+'-'+ platform
                if sampleID not in D[gene]:
                    D[gene][sampleID] = [float(fields[i]), float(fields[i+1])]
                    if sampleID not in Sample:
                        Sample.append(sampleID)
                else:
                    print(sampleID)
        inFile.close()

ouFile = open('TCGA_RNA-Seq_GeneCounts', 'w')
ouFile.write('Gene' + '\t' + '\t'.join(Sample) + '\n')
for g in Gene:
    L = []
    for s in Sample:
        L.append(str(int(D[g][s][0])))
    ouFile.write(g + '\t' + '\t'.join(L) + '\n')
ouFile.close()
ouFile = open('TCGA_RNA-Seq_Genes', 'w')
for g in Gene:
    ouFile.write(g + '\n')
ouFile.close()
ouFile = open('TCGA_RNA-Seq_Samples', 'w')
for s in Sample:
    ouFile.write(s + '\n')
ouFile.close()
        
