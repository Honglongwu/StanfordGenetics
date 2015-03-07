def CancerSampleNumber():
    Sample = {}
    inFile = open('TCGA-Cosmic')
    D = {}
    for line in inFile:
        line = line.strip()
        fields = line.split('\t')
        cancer_type = fields[7]
        sample_id = fields[4]
        D.setdefault(cancer_type, [])
        D[cancer_type].append(sample_id)
    inFile.close()
    for k in D:
        Sample[k] = len(set(D[k]))
    return Sample
CancerSample = CancerSampleNumber()
#inFile = open('NGLY1-TCGA-Cosmic')
inFile = open('TCGA-Cosmic')
ouFile = open('TCGA-Cosmic-Table', 'w')
ouFile2 = open('TCGA-Cosmic-Count', 'w')
#ouFile = open('NGLY1-TCGA-Cosmic-Table', 'w')
#ouFile2 = open('NGLY1-TCGA-Cosmic-Count', 'w')
D = {}
for line in inFile:
    line = line.strip()
    fields = line.split('\t')
    gene = fields[0]
    transcript = fields[1]
    cancer = fields[7]
    sample = fields[4]
    mutation = fields[15]
    homhet = fields[16]
    if not homhet:
        homhet = 'unknown'
    ouFile.write('\t'.join([gene, transcript, cancer, mutation, homhet, sample]) + '\n')


    D.setdefault(gene, {})
    D[gene].setdefault(transcript, {})
    D[gene][transcript].setdefault(cancer, {})
    D[gene][transcript][cancer].setdefault(mutation, 0)
    if homhet == 'hom':
        #D[gene][transcript][cancer][mutation] += 2
        D[gene][transcript][cancer][mutation] += 1
    elif homhet == 'het':
        D[gene][transcript][cancer][mutation] += 1
    else:
        D[gene][transcript][cancer][mutation] += 1

inFile.close()
ouFile.close()

for gene in D:
    for transcript in D[gene]:
        for cancer in D[gene][transcript]:
            for mutation in D[gene][transcript][cancer]:
                ouFile2.write('\t'.join([gene, transcript, cancer, mutation, str(D[gene][transcript][cancer][mutation]), str(CancerSample[cancer])]) + '\n')
ouFile2.close()


