Frequency = 0.01

inFile = open('ESP-snp-indels')
ouFile = open('ESP-snp-indels-rare', 'w')
for line in inFile:
    line = line.strip()
    fields = line.split(' ')
    gene =  fields[12]
    transcript =  fields[13]
    mutation = fields[14]
    genotype = fields[10]
    fds = genotype.split('/')
    L = []
    for fd in fds:
        x = int(fd.split('=')[-1])
        L.append(x)
    count_alt = sum(L[0:-1])
    count_ref = L[-1]
    if count_alt/float(count_alt + count_ref) < Frequency or count_ref/float(count_alt + count_ref) < Frequency:
        ouFile.write(line + '\n')

inFile.close()
ouFile.close()
