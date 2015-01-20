D = {}
inFile = open('x')
for line in inFile:
    fields = line.split('\t')
    D[fields[4]] = 1
inFile.close()

print(len(D))
