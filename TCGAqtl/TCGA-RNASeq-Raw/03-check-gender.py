def gender(n):
    D = {}
    D2 = {}
    inFile = open('TCGA_RNA-Seq_GeneCounts')
    ouFile = open('TCGA_RNA-Seq_SampleGender', 'w')
    head = inFile.readline().strip().split('\t')
    for line in inFile:
        if line[0:6] == 'KDM5D|':
            line = line.strip()
            fields = line.split('\t')
            for i in  range(1,len(fields)):
                if int(fields[i]) > n:
                    D[head[i]] = 'MALE'
                else:
                    D[head[i]] = 'FEMALE'
                D2[head[i]] = int(fields[i])
                ouFile.write(head[i] + '\t' + fields[i] + '\n')
    inFile.close()
    
    D3 = {}
    inFile = open('TCGA_RNA-Seq_SampleAnnot')
    head = inFile.readline()
    for line in inFile:
        fields = line.split('\t')
        D3[fields[0]] = fields[4]
    inFile.close()
    L = [] 
    for k in D:
        if D3[k]!='None':
            if D[k] != D3[k]:
                L.append(k)
    return(len(L))
                #print(k + '\t' + str(D2[k]))
#for n in range(200):
#    print(str(n) + '\t' + str(gender(n)))
    
    
def gender2():
    D = {}
    D2 = {}
    inFile = open('TCGA_RNA-Seq_GeneCounts')
    ouFile = open('TCGA_RNA-Seq_SampleGender', 'w')
    head = inFile.readline().strip().split('\t')
    for line in inFile:
        if line[0:6] == 'KDM5D|':
            line = line.strip()
            fields = line.split('\t')
            for i in  range(1,len(fields)):
                if int(fields[i]) > 20:
                    D[head[i]] = 'MALE'
                else:
                    D[head[i]] = 'FEMALE'
                D2[head[i]] = int(fields[i])
                ouFile.write(head[i] + '\t' + fields[i] + '\n')
    inFile.close()
    
    D3 = {}
    inFile = open('TCGA_RNA-Seq_SampleAnnot')
    head = inFile.readline()
    for line in inFile:
        fields = line.split('\t')
        D3[fields[0]] = fields[4]
    inFile.close()
    L = []
    for k in D:
        if D3[k]=='None':
            L.append(D2[k])
    L.sort()
    for x in L:
        print(x)
gender2()
