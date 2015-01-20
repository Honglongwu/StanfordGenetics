D = {}
inFile = open('TCGA-Cosmic-Count')
ouFile = open('TCGA-Cosmic-Count-unified', 'w')
for line in inFile:
    line = line.strip()
    fields = line.split('\t')
    mutation = fields[3]
    mutation_type = ''
    if mutation == 'Deletion - Frameshift' or mutation == 'Insertion - Frameshift':
                mutation_type = 'Frameshift'
    elif mutation == 'Substitution - coding silent':
        mutation_type = 'Synonymous'
    elif mutation == 'Substitution - Missense':
        mutation_type = 'Missense'
    elif mutation == 'Substitution - Nonsense' or mutation == 'Nonstop extension':
        mutation_type = 'Nonsense'
    ouFile.write(line + '\t' + mutation_type + '\n') 
inFile.close()
ouFile.close()

