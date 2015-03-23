from scipy import stats
import matplotlib
matplotlib.use('Agg')
import pylab as pl
SIG = 0.01

def Annot():
    SampleAnnot = {}
    inFile = open('TCGA_RNA-Seq_SampleAnnot')
    head = inFile.readline().strip().split('\t')
    for line in inFile:
        line = line.strip()
        fields = line.split('\t')
        sample = fields[0]
        SampleAnnot.setdefault(sample, {})
        for i in range(1,len(fields)):
            SampleAnnot[sample][head[i]] = fields[i]
    inFile.close()
    return SampleAnnot
SampleAnnot = Annot()

def getSample(L,cancer, sample_type='tumor'):
    S1 = []
    S2 = []
    SampleAnnot = Annot()
    for s in SampleAnnot:
        if SampleAnnot[s]['sample'] in L:
            S1.append(s)
        if SampleAnnot[s]['cancer'] in cancer:
            if SampleAnnot[s]['sample_type'] == sample_type:
                if SampleAnnot[s]['sample'] not in L:
                    S2.append(s)
    return([S1,S2])

#UCEC= ['TCGA-D1-A103','TCGA-D1-A17Q','TCGA-B5-A11N','TCGA-B5-A0JY']
CSMD3=['TCGA-64-5815','TCGA-05-4405','TCGA-53-7624','TCGA-64-1676','TCGA-44-6779','TCGA-55-1596','TCGA-44-6144','TCGA-34-5927','TCGA-22-5491','TCGA-44-7662','TCGA-66-2783','TCGA-22-5472','TCGA-50-5946','TCGA-44-6775','TCGA-37-3789','TCGA-53-7813','TCGA-18-3414','TCGA-95-7039','TCGA-18-3409','TCGA-18-3409','TCGA-22-5471','TCGA-73-4670','TCGA-60-2698','TCGA-44-7670','TCGA-60-2711','TCGA-33-4582','TCGA-60-2711','TCGA-80-5611','TCGA-66-2795','TCGA-50-5946-02','TCGA-78-7161','TCGA-34-5231','TCGA-78-7163','TCGA-44-4112','TCGA-37-5819','TCGA-55-1595','TCGA-44-5643']
#S=getSample(UCEC,'UCEC')
S=getSample(CSMD3,['LUAD','LUSC'])

def boxplot(X,Y,gene,p):
    pl.figure()
    pl.boxplot([X,Y])
    pl.xticks([1,2],["mutated","control"])
    gs = gene.split('|')
    if gs[0] == '?':
        g = 'Unknown' + '_' + gs[1]
    else:
        g = gs[0] + '_' + gs[1]

    if gs[0] != '?':
        pl.title('%s (p-value = %s)'%(gs[0],p))
    pl.savefig('pdf/'+ g + '.sig.pdf')


def Sig(S1, S2):
    P = []
    inFile = open('TCGA_RNA-Seq_GeneCounts-Normalized')
    head = inFile.readline().strip().split('\t')
    for line in inFile:
        line = line.strip()
    #for i in range(500):
    #    line = inFile.readline().strip()
        X = []
        Y= []
        fields = line.split('\t')
        gene = fields[0]
        for i in range(1, len(fields)):
            if head[i] in S1:
                X.append(float(fields[i]))
            elif head[i] in S2:
                Y.append(float(fields[i]))
        #p = stats.ttest_ind(X, Y)[1]
        X.sort()
        Y.sort()
        p = stats.ranksums(X, Y)[1]
        P.append([gene, p, X, Y])
    inFile.close()
    P.sort(cmp = lambda x,y:cmp(x[1], y[1]))
    for m in P:
        if m[1] < SIG:
            print(m[0] + '\t' + str(m[1])+'\t' +'X=['+ ','.join([str(x) for x in m[2]])+']' + '\t' +'Y=['+ ','.join([str(x) for x in m[3]])+']')
            boxplot(m[2], m[3], m[0], m[1])
        else:
            print(m[0] + '\t' + str(m[1]))
Sig(S[0], S[1])






