def classification(cancer):
    inFile = open('TCGA_RNA-Seq_GeneCounts')
    head = inFile.readline().strip().split('\t')
    ouFile = open('TCGA_RNA-Seq_GeneCounts-' + cancer,'w')
    H = [head[0]]
    for i in range(len(head)):
        if cancer in head[i] and head[i].split('-')[3][0] == '0':
            H.append(head[i])
    ouFile.write('\t'.join(H) + '\n')
    for line in inFile:
        line = line.strip()
        fields = line.split('\t')
        L = [fields[0]]
        for i in range(len(fields)):
            if cancer in head[i]:
                if head[i].split('-')[3][0] == '0':
                    L.append(fields[i])
        ouFile.write('\t'.join(L) + '\n')
    inFile.close()
    ouFile.close()

classification('UCEC')
