inFile = open('TCGA-ESP-count')
ouFile = open('TCGA-ESP-count-enriched', 'w')
for line in inFile:
    line = line.strip()
    fields = line.split('\t')
    freq1 = float(fields[4])/float(fields[5])
    freq2 = float(fields[10])/float(fields[11])
    if freq1 > freq2:
        ouFile.write(line + '\n')

inFile.close()
ouFile.close()
