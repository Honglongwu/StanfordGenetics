import os

def clinicalData():
    D = {}
    DIR = '/srv/gsfs0/projects/steinmetz/hansun/TCGAqtl/TCGA-Clinical'
    Fs = os.listdir(DIR)
    for F in Fs:
        if F[-4:] == '.txt':
            inFile = open(DIR + os.path.sep + F)
            head1 = inFile.readline().strip('\n').split('\t')
            head2 = inFile.readline().strip('\n').split('\t')
            head3 = inFile.readline().strip('\n').split('\t')
            for line in inFile:
                line = line.strip('\n')
                if line:
                    fields = line.split('\t')
                    sample = fields[head1.index('bcr_patient_barcode')]
                    D.setdefault(sample, {})
                    for i in range(len(fields)):
                        if fields[i] == '[Not Available]' or fields[i] == '[Not Applicable]':
                            D[sample][head1[i]] = 'None'
                        else:
                            D[sample][head1[i]] = fields[i]
            inFile.close()
    for k in D:
        if 'age_at_diagnosis' in D[k]:
            age = D[k]['age_at_diagnosis']
        elif 'age_at_initial_pathologic_diagnosis' in D[k]:
            age = D[k]['age_at_initial_pathologic_diagnosis']
        else:
            age='None'
        D[k]['age'] = age
        if 'tumor_tissue_site' in D[k]:
            tissue = D[k]['tumor_tissue_site']
        elif 'tumor_tissue_site_other' in D[k]:
            tissue = D[k]['tumor_tissue_site_other']
        else:
            tissue='None'
        D[k]['tissue'] = tissue

        if 'tumor_status' not in D[k]:
            tumor_status='None'
            D[k]['tumor_status'] = tumor_status

        #print('\t'.join([D[k]['bcr_patient_barcode'],D[k]['gender'],D[k]['race'],age,D[k]['tumor_status'], tissue]))
    return D

def gender():
    D = {}
    inFile = open('TCGA_RNA-Seq_GeneCounts')
    head = inFile.readline().strip().split('\t')
    for line in inFile:
        if line[0:6] == 'KDM5D|':
            line = line.strip()
            fields = line.split('\t')
            for i in  range(1,len(fields)):
                sn = '-'.join(head[i].split('-')[0:3])
                if int(fields[i]) > 20:
                    D[sn] = 'MALE'
                else:
                    D[sn] = 'FEMALE'
    inFile.close()
    return D


D = {}
Samples = []
CL = clinicalData()
inFile = open('TCGA_RNA-Seq_Samples')
for line in inFile:
    line = line.strip()
    Samples.append(line)
inFile.close()

for s in Samples:
    sf = s.split('-')
    D.setdefault(s, {})
    D[s]['sample'] = '-'.join(sf[0:3])
    D[s]['cancer'] = sf[-2]
    D[s]['platform'] = sf[-1]
    #if D[s]['sample'] not in CL:
    #    print(s)
G = gender()
ouFile = open('TCGA_RNA-Seq_SampleAnnot', 'w')
ouFile.write('\t'.join(['sampleID', 'sample', 'cancer', 'platform','gender','race','age','tumor_status','tissue']) + '\n')
for s in Samples:
    cl = []

    sn = D[s]['sample']
    if sn in CL:
        cl = [CL[sn]['gender'], CL[sn]['race'], CL[sn]['age'], CL[sn]['tumor_status'], CL[sn]['tissue']]
    else:
        cl = ['None','None','None','None','None']
    
    if cl[0] == 'None':    
        cl[0] = G[sn]
        print('\t'.join([s, D[s]['sample'], D[s]['cancer'], D[s]['platform']]+cl))
    ouFile.write('\t'.join([s, D[s]['sample'], D[s]['cancer'], D[s]['platform']]+cl) + '\n')
ouFile.close()
