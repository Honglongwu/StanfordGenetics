def genotypeFormat(inF, c):
    """
    NA or N/A: -1
    WA: 0
    SA: 1
    WE: 2

    g012 to hdf5:
    limix_converter --outfile=Yeast-Genotype.hdf5 --g012=Yeast-Genotype.012

    following strains do not have phenotype:
    Mat_a_37
    Mat_a_38
    Mat_a_39
    Mat_a_40
    Mat_alpha_37
    Mat_alpha_38
    Mat_alpha_39
    Mat_alpha_40

    
    """
    Ex = ['Mat_a_37','Mat_a_38','Mat_a_39','Mat_a_40','Mat_alpha_37','Mat_alpha_38','Mat_alpha_39','Mat_alpha_40']
    inFile = open(inF)
    L = []
    Sample = []
    head = inFile.readline().rstrip()
    Position = head.split('\t')[1:]
    for line in inFile:
        line = line.rstrip()
        fields = line.split('\t')
        xL = []
        if fields[0] not in Ex:
            Sample.append(fields[0])
            for item in fields[1:]:
                if item == 'N/A':
                    xL.append(str(float('NaN')))
                elif item == c:
                    xL.append('1')
                else:
                    xL.append('0')
            L.append(xL)
    inFile.close()
    
    ouFile = open(inF.split('.')[0] +'-'+c+ '.012', 'w')
    n = 0
    for item in L:
        ouFile.write(str(n) + '\t' + '\t'.join(item) + '\n')
        n += 1
    ouFile.close()
    ouFile = open(inF.split('.')[0] +'-'+c+ '.012.indv', 'w')
    ouFile.write('\n'.join(Sample) + '\n')
    ouFile.close()
    ouFile = open(inF.split('.')[0] +'-'+c+ '.012.pos', 'w')
    ouFile.write('\n'.join(Position) + '\n')
    ouFile.close()

genotypeFormat('Yeast-Genotype-noMissing.tab', 'WA')
genotypeFormat('Yeast-Genotype-noMissing.tab', 'SA')
genotypeFormat('Yeast-Genotype-noMissing.tab', 'NA')
genotypeFormat('Yeast-Genotype-noMissing.tab', 'WE')
