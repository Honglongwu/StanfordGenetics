inFile = open('ESP-snp-indels')
#inFile = open('NGLY1-ESP-snp-indels')
#ouFile = open('NGLY1-ESP-snp-indels-Count', 'w')
ouFile = open('ESP-snp-indels-Count', 'w')

D = {}
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
    D.setdefault(gene, {})
    D[gene].setdefault(transcript, {})
    D[gene][transcript].setdefault(mutation,[0])
    #D[gene][transcript].setdefault(mutation,[0,0])
    D[gene][transcript][mutation][0] += count_alt
    #D[gene][transcript][mutation][1] += count_ref
    #D[gene][transcript][mutation].append(count_ref)
    D[gene][transcript][mutation].append(count_ref + count_alt)
inFile.close()

for gene in D:
    for transcript in D[gene]:
        for mutation in D[gene][transcript]:
            #ouFile.write('\t'.join([gene, transcript, mutation, str(D[gene][transcript][mutation][0]), str(max(D[gene][transcript][mutation][1:])), str(min(D[gene][transcript][mutation][1:])) ]) + '\n')
            ouFile.write('\t'.join([gene, transcript, mutation, str(D[gene][transcript][mutation][0]), str(max(D[gene][transcript][mutation][1:]))]) + '\n')
            #ouFile.write('\t'.join([gene, transcript, mutation, str(D[gene][transcript][mutation][0]), str(D[gene][transcript][mutation][1])]) + '\n')
ouFile.close()
