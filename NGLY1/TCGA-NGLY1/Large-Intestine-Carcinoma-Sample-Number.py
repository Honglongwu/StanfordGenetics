D = {}
D2 = {}
D3 = {}

inFile = open('READ__unc.edu__illuminaga_rnaseqv2__rsem.genes.results__Jul-08-2014.txt')
head = inFile.readline().strip().split('\t')
for x in head[2:]:
    D[x]=1
    D2[x]=1
inFile.close()

inFile = open('READ__unc.edu__illuminahiseq_rnaseqv2__rsem.genes.results__Jul-08-2014.txt')
head = inFile.readline().strip().split('\t')
for x in head[2:]:
    D[x]=1
    D2[x]=1
inFile.close()

inFile = open('COAD__unc.edu__illuminaga_rnaseqv2__rsem.genes.results__Jul-08-2014.txt')
head = inFile.readline().strip().split('\t')
for x in head[2:]:
    D[x]=1
    D3[x]=1
inFile.close()

inFile = open('COAD__unc.edu__illuminahiseq_rnaseqv2__rsem.genes.results__Jul-08-2014.txt')
head = inFile.readline().strip().split('\t')
for x in head[2:]:
    D[x]=1
    D3[x]=1
inFile.close()

print(len(D))
print(len(D2))
print(len(D3))
