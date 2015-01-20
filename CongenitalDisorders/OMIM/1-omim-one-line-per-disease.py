inFile = open('omim.txt')
ouFile = open('omim-one-line-per-disease', 'w')
L = []
headi = -1
for line in inFile:
    line = line.strip()
    if line.find('*RECORD*') != -1:
        headi += 1
        L.append([line])
    else:
        L[headi].append(line)
inFile.close()
for item in L:
    ouFile.write('\t'.join(item)+'\n')
ouFile.close()
