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

READ_COAD= ['TCGA-AG-A002','TCGA-EI-6917','TCGA-A6-6141']
S=getSample(READ_COAD,['READ','COAD'])

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
    pl.savefig('pdf/'+ g + '.READ_COAD.sig.pdf')


def Sig(S1, S2):
    P = []
    ouFile = open('READ_COAD-Sig.txt', 'w')
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
            ouFile.write(m[0] + '\t' + str(m[1])+'\t' +'X=['+ ','.join([str(x) for x in m[2]])+']' + '\t' +'Y=['+ ','.join([str(x) for x in m[3]])+']'+ '\n')
            boxplot(m[2], m[3], m[0], m[1])
        else:
            ouFile.write(m[0] + '\t' + str(m[1]) + '\n')

    ouFile.close()
Sig(S[0], S[1])







