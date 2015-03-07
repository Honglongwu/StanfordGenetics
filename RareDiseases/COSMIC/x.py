inFile = open('CosmicCompleteExport.tsv')
D = {}
for line in inFile:
    fields = line.split('\t')
    D[fields[15]] =1 
inFile.close()

for  k in D:
    print(k)
