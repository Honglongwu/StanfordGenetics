D = {}
D2 = {}
inFile = open("nationwidechildrens.org_clinical_patient_ucec.txt")
head = inFile.readline()
head = inFile.readline()
head = inFile.readline()

for line in inFile:
    line = line.strip()
    fields = line.split("\t")
    if len(fields) > 1:
        sample = fields[0]
        D[sample] = [fields[6], fields[17], fields[21], fields[31]]
inFile.close()

inFile = open("UCEC-samples")
ouFile = open("UCEC-samples-phenotype", 'w')
for line in inFile:
    fields = line.split()
    fds = fields[0].split('-')
    k = '-'.join(fds[0:3])
    if k in D:
        ouFile.write(fields[0] + '\t' + '\t'.join(D[k])+'\n')
    else:
        ouFile.write(fields[0] + '\t' + '\t'.join(['NA','NA','NA','NA'])+'\n')
inFile.close()
ouFile.close()

