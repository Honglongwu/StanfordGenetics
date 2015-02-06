NM2ENST = {}
inFile = open('/srv/gsfs0/projects/steinmetz/hansun/Data/gene2ensembl-human')
for line in inFile:
    line = line.strip()
    fields = line.split('\t')
    NM = fields[3]
    ENST = fields[4]
    NM2ENST[NM] = ENST
    #NM2ENST[NM.split('.')[0]] = ENST
inFile.close()
D = {}
inFile = open('ESP-snp-indels-Count')
ouFile = open('ESP-snp-indels-Count-unified', 'w')
for line in inFile:
    line = line.strip()
    fields = line.split('\t')
    mutation = fields[2]
    mutation_type = ''
    enst = ''
    nm =  fields[1]
    if mutation == 'frameshift':
        mutation_type = 'Frameshift'
    elif mutation == 'coding-synonymous':
        mutation_type = 'Synonymous'
    elif mutation == 'missense':
        mutation_type = 'Missense'
    elif mutation == 'stop-lost' or mutation == 'stop-gained':
        mutation_type = 'Nonsense'
    if nm in NM2ENST:
        enst = NM2ENST[nm]
    #elif nm.split('.')[0] in NM2ENST:
    #    enst = NM2ENST[nm.split('.')[0]]
    #    print(line)
    ouFile.write(line + '\t' + mutation_type + '\t' + enst +'\n') 
inFile.close()
ouFile.close()

