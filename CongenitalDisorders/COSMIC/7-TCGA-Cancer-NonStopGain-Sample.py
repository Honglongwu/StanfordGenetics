inFile = open('TCGA-Sample-NonStopGain-Cancer')
ouFile = open('TCGA-Cancer-NonStopGain-Sample-Frequency', 'w')
D = {}
for line in inFile:
    line = line.strip()
    fields = line.split('\t')
    D[fields[0]] = fields[1]
inFile.close()

inFile = open('TCGA-Cosmic-Sample-Number')
for line in inFile:
    line = line.strip()
    fields = line.split('\t')
    cancer = fields[0]
    number = fields[1]
    frequency = '%.2f' % (float(D.get(cancer, 0))/float(number))
    ouFile.write(cancer + '\t' + str(D.get(cancer, 0)) + '\t' + number + '\t' +frequency +'\n')
inFile.close()
ouFile.close()
