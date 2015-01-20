D = {}
inFile = open('ESP-snp-indels-Count')
ouFile = open('ESP-snp-indels-Count-unified', 'w')
for line in inFile:
    line = line.strip()
    fields = line.split('\t')
    mutation = fields[2]
    mutation_type = ''
    if mutation == 'frameshift':
        mutation_type = 'Frameshift'
    elif mutation == 'coding-synonymous':
        mutation_type = 'Synonymous'
    elif mutation == 'missense':
        mutation_type = 'Missense'
    elif mutation == 'stop-lost' or mutation == 'stop-gained':
        mutation_type = 'Nonsense'
    ouFile.write(line + '\t' + mutation_type + '\n') 
inFile.close()
ouFile.close()

