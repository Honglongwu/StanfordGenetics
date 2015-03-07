inFile = open('ESP-snp-indels')
for line in inFile:
    line = line.strip()
    fields = line.split(' ')
    p = fields[25]
    if p == 'PASS':
        pass
    else:
        print(line)

inFile.close()
