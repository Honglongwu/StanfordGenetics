inFile = open('omim.txt')
ouFile = open('omim-one-line-per-disease-withoutREFENCE', 'w')
L = []
headi = -1
RFflag = 0
for line in inFile:
    line = line.strip()
    if line.find('*RECORD*') != -1:
        headi += 1
        L.append([line])
    else:
        if line.find('*FIELD* RF') != -1:
            RFflag = 1
        elif line.find('*') == 0:
            RFflag = 0
        if RFflag == 0:
            L[headi].append(line)
inFile.close()
for item in L:
    ouFile.write('\t'.join(item)+'\n')
ouFile.close()
