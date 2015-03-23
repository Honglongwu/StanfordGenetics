def phenotypeFormat(inF):
    inFile = open(inF)
    head = inFile.readline().rstrip('\n')
    D = {}
    ouFile2 = open(inF.split('.txt')[0] + '.strain.txt', 'w')
    PS = []
    for line in inFile:
        line = line.rstrip('\n')
        fields = line.split('\t')
        environment = fields[0]
        component = fields[1]
        strain = fields[2]
        PS.append(strain)
        lsc_mean = fields[9]
        lsc_sd = fields[10]
        D.setdefault(strain, {})
        D[strain].setdefault(component, {})
        D[strain][component].setdefault(environment, {})
        D[strain][component][environment]=[lsc_mean, lsc_sd]
    inFile.close()
    ouFile2.write('\n'.join(set(PS)) + '\n')
    ouFile2.close()
    
    E = []
    C = []
    for strain in D:
        for component in D[strain]:
            for environment in D[strain][component]:
                E.append(environment)
            break
        break
    
    for strain in D:
        for component in D[strain]:
            C.append(component)
    C = set(C)
    
    ouFile = open(inF.split('.txt')[0]+'.formated.txt', 'w')
    #### output head
    head = []
    for c in C:
        for e in E:
            head.append('%s-%s-mean\t%s-%s-sd'%(c,e,c,e))
    ouFile.write('Strain\t' + '\t'.join(head) + '\n')
    
    S = strainList()
    for strain in S:
        strainN = strainName(strain)
        if strainN:
            if strainN in D:
                ouFile.write(strain + '\t')
                CE = []
                for component in D[strainN]:
                    for environment in E:
                        #ouFile.write('\t'.join(D[strain][component][environment])+'\t')
                        CE += D[strainN][component][environment]
                ouFile.write('\t'.join(CE) + '\n')
            else:
                #print(strain)
                pass
        else:
            #print(strain)
            pass

    ouFile.close()

def strainList():
    L = []
    inFile = open('Yeast-Genotype-noMissing-SA.012.indv')
    for line in inFile:
        line = line.strip()
        L.append(line)
    inFile.close()
    return L

def checkStrain():
    L1 = []
    L2 = []
    inFile = open('Yeast-Genotype-noMissing-SA.012.indv')
    for line in inFile:
        line = line.strip()
        L1.append(line)
    inFile.close()

    inFile = open('Yeast-Phenotype.txt')
    head = inFile.readline()
    for line in inFile:
        line = line.strip()
        fields = line.split('\t')
        L2.append(fields[2])
    inFile.close()
    L2 = set(L2)
    print('Warning - strain has genotype, but not phenotype:')
    for x in L1:
        if strainName(x) not in L2:
            print(x)
    print('Warning - strain has phenotype, but not genotype:')
    for x in L2:
        if strainName2(x) not in L1:
            print(x)


def strainName2(strain):
    name = ''
    if strain[0:2] == '20':
        name = 'Mat_alpha_' + strain[2:] 
    elif strain[0] == '2':
        name = 'Mat_alpha_' + strain[1:] 
    elif strain[0:2] == '30':
        name = 'Mat_a_' + strain[2:] 
    elif strain[0] == '3':
        name = 'Mat_a_' + strain[1:] 
    return name
def strainName(strain):
    name = ''
    if strain[0:10] == 'Mat_alpha_':
        name = strain[10:] 
        if len(name) == 1:
            name = '20' + name 
        else:
            name = '2' + name 
    elif strain[0:6] == 'Mat_a_':
        name = strain[6:] 
        if len(name) == 1:
            name = '30' + name 
        else:
            name = '3' + name 
    return name


phenotypeFormat('Yeast-Phenotype.txt')
checkStrain()

