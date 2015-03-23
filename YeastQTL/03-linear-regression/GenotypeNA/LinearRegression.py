import matplotlib
matplotlib.use('Agg')
import limix.io.genotype_reader as gr
import limix.io.phenotype_reader as phr
import limix.io.data as data
import scipy as sp
import pylab as pl
import pandas as pd
import limix.io.data_util as data_util
from limix.utils.plot import *
import limix.utils.preprocess as preprocess
import limix.modules.qtl as qtl
import limix.stats.fdr as fdr

geno_reader = gr.genotype_reader_tables('../Yeast-Genotype-noMissing-NA.hdf5')
pheno_reader = phr.pheno_reader_tables('../Yeast-Phenotype.hdf5')
dataset = data.QTLData(geno_reader=geno_reader,pheno_reader=pheno_reader)
geno = dataset.getGenotypes()
position = dataset.getPos()
pos,chromBounds = data_util.estCumPos(position=position,offset=0)
P_max = len(dataset.phenotype_ID)
#P_max = 6
phenotype_ID = dataset.phenotype_ID[0:P_max:2]
phenotype_vals, sample_idx = dataset.getPhenotypes(phenotype_ID)
N = geno.shape[0]
S = geno.shape[1]
P = phenotype_vals.shape[1]

def phenoPlot(phenotype_ID, phenotype_vals):
    for ip, p_ID in enumerate(phenotype_ID):
        pl.figure()
        plot_normal(phenotype_vals.values[:,ip],alpha=0.8)
        pl.title("%s" % p_ID)
        pl.savefig('Phenotype-%s.pdf'%p_ID)
        pl.close('all')

phenoPlot(phenotype_ID, phenotype_vals)

phenotype_vals_boxcox, maxlog = preprocess.boxcox(phenotype_vals.values)
phenotype_vals_ranks = preprocess.rankStandardizeNormal(phenotype_vals.values)

def phenoCmpPlot(phenotype_ID, phenotype_vals):
    for ip, p_ID in enumerate(phenotype_ID):
        pl.figure(figsize=[10,3])#create the figure
        
        plt = pl.subplot(1,3,1)#the untransformed phenotypes
        #histogram of the untransformed phenotypes
        plot_normal(phenotype_vals.values[:,ip],alpha=0.8,figure=plt)
        pl.title("untransformed")
        
        plt = pl.subplot(1,3,2)#the untransformed phenotypes
        #histogram of the untransformed phenotypes
        plot_normal(phenotype_vals_boxcox[:,ip],alpha=0.8,figure=plt)
        #pl.title("%s, Box-Cox" % p_ID)
        pl.title("Box-Cox")
        
        plt = pl.subplot(1,3,3)#the rank transformed phenotypes
        #histogram of the untransformed phenotypes
        plot_normal(phenotype_vals_ranks[:,ip],alpha=0.8,figure=plt)
        pl.title("ranks")
        pl.savefig('Phenotype-transformation-%s.pdf'%p_ID)
        pl.close('all')

phenoCmpPlot(phenotype_ID, phenotype_vals)

SIG = 0.000001
#SIG = 0.0001
def significant(pvalues_lm,ouF):
    SigPhe = set()
    ALL = []
    ouFile = open(ouF,'w')
    ouFile2 = open(ouF+'-ALL','w')
    for i in range(pvalues_lm.shape[1]):
        n = -1
        for x in pvalues_lm.ix[:,i]:
            n += 1
            ALL.append([n, pos['chrom'][n], pos['pos'][n], i*2, pvalues_lm.columns[i], x])
            if x < SIG:
                ouFile.write("%s\t%s_%s\t%s\t%s\t%s"%(n,pos['chrom'][n], pos['pos'][n],i*2,pvalues_lm.columns[i],x) + '\n')
                SigPhe.add(pvalues_lm.columns[i])
    ouFile.close()
    ALL.sort(cmp = lambda x,y:cmp(x[-1],y[-1]))
    for x in ALL:
        ouFile2.write("%s\t%s_%s\t%s\t%s\t%s"%(x[0],x[1],x[2],x[3],x[4],x[5]) + '\n')
    ouFile2.close()
    return SigPhe

def manhattonPlot(phenotype_ID, pvalues_lm, ouFprefix):
    for ip, p_ID in enumerate(phenotype_ID):
        pl.figure(figsize=[12,4])
        plot_manhattan(posCum=pos['pos_cum'],pv=pvalues_lm[p_ID].values,chromBounds=chromBounds,thr_plotting=0.05)
        pl.title(p_ID)
        pl.savefig(ouFprefix + '.' + p_ID + '.pdf')
        pl.close('all')


lm = qtl.test_lm(snps=geno[sample_idx],pheno=phenotype_vals.values)
pvalues_lm = pd.DataFrame(data=lm.pvalues.T,index=dataset.geno_ID,columns=phenotype_ID)
SigPhe = significant(pvalues_lm, 'QTL-lm-significant-untransformed')
manhattonPlot(SigPhe, pvalues_lm, 'QTL-lm-significant-untransformed')


lm_boxcox = qtl.test_lm(snps=geno[sample_idx],pheno=phenotype_vals_boxcox)
pvalues_lm_boxcox = pd.DataFrame(data=lm_boxcox.pvalues.T,index=dataset.geno_ID,columns=phenotype_ID)
SigPhe = significant(pvalues_lm_boxcox, 'QTL-lm-significant-boxcox')
manhattonPlot(SigPhe, pvalues_lm_boxcox, 'QTL-lm-significant-boxcox')


lm_ranks = qtl.test_lm(snps=geno[sample_idx],pheno=phenotype_vals_ranks)
pvalues_lm_ranks = pd.DataFrame(data=lm_ranks.pvalues.T,index=dataset.geno_ID,columns=phenotype_ID)
SigPhe = significant(pvalues_lm_ranks,'QTL-lm-significant-rank')
manhattonPlot(SigPhe, pvalues_lm_ranks, 'QTL-lm-significant-rank')

def transformCmp(pvalues_lm, pvalues_lm_boxcox, pvalues_lm_ranks, ouFprefix):
    for ip, p_ID in enumerate(dataset.phenotype_ID[0:P_max:2]):
        #pl.figure(figsize=[15,3])#create the figure
        pl.figure(figsize=[15,5])#create the figure
        pl.subplot(1,3,1)
        pl.plot(-sp.log10(pvalues_lm[p_ID].values),-sp.log10(pvalues_lm_ranks[p_ID].values),'.')
        pl.xlabel('$-log_{10}(pv)$ lm')
        pl.ylabel('$-log_{10}(pv)$ lm ranks')
        pl.title(p_ID)
        max_range = max(pl.xlim()[1],pl.ylim()[1])
        pl.plot([0,max_range],[0,max_range],'k--')
        pl.subplot(1,3,2)
        pl.plot(-sp.log10(pvalues_lm[p_ID].values),-sp.log10(pvalues_lm_boxcox[p_ID].values),'.')
        pl.xlabel('$-log_{10}(pv)$ lm')
        pl.ylabel('$-log_{10}(pv)$ lm boxcox')
        max_range = max(pl.xlim()[1],pl.ylim()[1])
        pl.plot([0,max_range],[0,max_range],'k--')
        pl.title(p_ID)
        pl.subplot(1,3,3)
        pl.plot(-sp.log10(pvalues_lm_boxcox[p_ID].values),-sp.log10(pvalues_lm_ranks[p_ID].values),'.')
        pl.xlabel('$-log_{10}(pv)$ lm boxcox')
        pl.ylabel('$-log_{10}(pv)$ lm ranks')
        pl.title(p_ID)
    
        max_range = max(pl.xlim()[1],pl.ylim()[1])
        pl.plot([0,max_range],[0,max_range],'k--')
        pl.savefig(ouFprefix + '.' + p_ID + '.pdf')
        pl.close('all')

transformCmp(pvalues_lm, pvalues_lm_boxcox, pvalues_lm_ranks, 'QTL-lm-transformation-comparison')

def phenoGeno(lm, ouFprefix):
    for ip, p_ID in enumerate(dataset.phenotype_ID[0:P_max:2]):
        pl.figure(figsize=[15,5])#create the figure

        plt = pl.subplot(1,3,1)#the untransformed phenotypes
        #find maximum squared beta value
        pheno_vals, s_idx = dataset.getPhenotypes([p_ID])
        
        imax = lm.pvalues[ip].argmin()
        i_0 = geno[s_idx,imax]==0
        
        #plot SNP vs. phenotype for max beta
        
        pl.plot(geno[s_idx,imax]+0.05*np.random.randn(geno[s_idx,imax].shape[0]),pheno_vals.values,'.',alpha=0.5)
        
        pl.xlabel("Genotyp")
        pl.ylabel("Phenotype")
        pl.xlim([-0.5,1.5])
        pl.title("%s" % p_ID)
        
        plt = pl.subplot(1,3,2)#the Box-Cox transformed phenotypes
        pl.plot(geno[s_idx,imax]+0.05*np.random.randn(geno[s_idx,imax].shape[0]),phenotype_vals_boxcox[s_idx[sample_idx],ip],'.',alpha=0.5)
        
        pl.xlabel("Genotype")
        pl.ylabel("Phenotype")
        pl.xlim([-0.5,1.5])
        pl.title("%s, Box-Cox" % p_ID)
        
        plt = pl.subplot(1,3,3)#the rank transformed phenotypes
        pl.plot(geno[s_idx,imax]+0.05*np.random.randn(geno[s_idx,imax].shape[0]),phenotype_vals_ranks[s_idx[sample_idx],ip],'.',alpha=0.5)
        
        #pl.plot([0,1],[pheno_vals.values[i_0].mean(),pheno_vals.values[~i_0].mean()])
        pl.xlabel("Genotype")
        pl.ylabel("Phenotype")
        pl.xlim([-0.5,1.5])
        pl.title("%s, ranks" % p_ID)

        pl.savefig(ouFprefix + '.' + p_ID + '.pdf')
        pl.close('all')
    
phenoGeno(lm, 'QTL-lm-genotype-phenotype')

def pvaluesHist(pvalues_lm, pvalues_lm_boxcox, pvalues_lm_ranks, ouFprefix):
    for ip, p_ID in enumerate(dataset.phenotype_ID[0:P_max:2]):
        pl.figure(figsize=[15,5])

        plt = pl.subplot(1,3,1)
        pl.hist(pvalues_lm[p_ID].values,20,normed=True)
        pl.plot([0,1],[1,1],"r")
        pl.title("%s" % p_ID)
        pl.xlabel("P-value")
        pl.ylabel("Frequency")
        
        plt = pl.subplot(1,3,2)
        pl.hist(pvalues_lm_boxcox[p_ID].values,20,normed=True)
        pl.plot([0,1],[1,1],"r")
        pl.title("%s, Box-Cox" % p_ID)
        pl.xlabel("P-value")
        
        plt = pl.subplot(1,3,3)
        pl.hist(pvalues_lm_ranks[p_ID].values,20,normed=True)
        pl.plot([0,1],[1,1],"r")
        pl.title("%s, ranks" % p_ID)
        pl.xlabel("P-value")

        pl.savefig(ouFprefix + '.' + p_ID + '.pdf')
        pl.close('all')
pvaluesHist(pvalues_lm, pvalues_lm_boxcox, pvalues_lm_ranks, 'QTL-lm-pvalues-histgram')

def qqPlot(pvalues_lm, pvalues_lm_boxcox, pvalues_lm_ranks,ouFprefix):
    for ip, p_ID in enumerate(dataset.phenotype_ID[0:P_max:2]):
        pl.figure(figsize=[15,5])
        
        plt = pl.subplot(1,3,1)
        qqplot(pvalues_lm[p_ID].values)
        pl.title("%s" % p_ID)
        
        plt = pl.subplot(1,3,2)
        qqplot(pvalues_lm_boxcox[p_ID].values)
        pl.title("%s, Box-Cox" % p_ID)
        
        plt = pl.subplot(1,3,3)
        qqplot(pvalues_lm_ranks[p_ID].values)
        pl.title("%s, ranks" % p_ID)

        pl.savefig(ouFprefix + '.' + p_ID + '.pdf')
        pl.close('all')

qqPlot(pvalues_lm, pvalues_lm_boxcox, pvalues_lm_ranks, 'QTL-lm-pvalues-qqplot')

def permutationPlot(ouFprefix):
    phenotype_vals_perm = phenotype_vals.copy()

    for ip, p_ID in enumerate(dataset.phenotype_ID[0:P_max:2]):
        perm = np.random.permutation(phenotype_vals[p_ID].values)
        phenotype_vals_perm[p_ID] = perm
      
        #run linear regression on each SNP on permutated and non-permuted data
        lm = qtl.test_lm(snps=geno[sample_idx],pheno=phenotype_vals.values)
        lm_perm = qtl.test_lm(snps=geno[sample_idx],pheno=phenotype_vals_perm.values)
        
        #plot
        fig = pl.figure(figsize=[12,5])
        pl.subplot(1,2,1)
        pl.hist(lm_perm.test_statistics.flatten(),normed=True,color='g',alpha=0.5,bins=50,range=(0,20),label='perm')
        pl.hist(lm.test_statistics[1],normed=True,color='b',alpha=0.5,bins=50,range=(0,20),label='observed')
        pl.xlabel('test statistics')
        pl.ylabel('density')
        pl.title(p_ID)
    
        pl.legend()
        pl.subplot(1,2,2)
        p=pl.hist(lm_perm.pvalues.flatten(),normed=True,color='g',alpha=0.5,bins=50,range=(0,1),label='perm')
        p=pl.hist(lm.pvalues.flatten(),normed=True,color='b',alpha=0.5,bins=50,range=(0,1),label='observed')
        pl.title(p_ID)
        pl.xlabel('p-value')
        pl.ylabel('density')
        pl.legend()
        pl.savefig(ouFprefix + '.' + p_ID + '.pdf')
        pl.close('all')

permutationPlot('QTL-lm-permutation')

def fdrPlot(ouFprefix):
    phenotype_vals_perm = phenotype_vals.copy()

    for ip, p_ID in enumerate(dataset.phenotype_ID[0:P_max:2]):
        perm = np.random.permutation(phenotype_vals[p_ID].values)
        phenotype_vals_perm[p_ID] = perm
      
        #run linear regression on each SNP on permutated and non-permuted data
        lm = qtl.test_lm(snps=geno[sample_idx],pheno=phenotype_vals.values)
        lm_perm = qtl.test_lm(snps=geno[sample_idx],pheno=phenotype_vals_perm.values)
 
        fig = pl.figure(ip, figsize=[15,5])
        
        #1. on the real tests
        pl.subplot(1,2,1)
        qv, qv_info = fdr.qvalues(lm.pvalues[1])
        pi0 = qv_info['pi0']
        #plot distribution of p-values and fitted \pi parameter:
        fig = pl.figure(ip,figsize=[8,5])
        pl.hist(lm.pvalues[1],bins=100,normed=True,color='b',label='observed',alpha=0.5)
        pl.axhline(pi0,color='g',linewidth=2,label='$\pi$')
        pl.axhline(1,color='r',linestyle='--',linewidth=2,label='null')
        pl.legend()
        pl.xlabel('p-value')
        pl.ylabel('density')
        pl.ylim((0,10))
        
        #2. on permuted p-values
        pl.subplot(1,2,2)
        qv_perm, qv_perm_info = fdr.qvalues(lm_perm.pvalues[1])
        pi0 = qv_perm_info['pi0']
        #plot distribution of p-values and fitted \pi parameter:
        pl.hist(lm_perm.pvalues[1],bins=100,normed=True,color='b',label='observed',alpha=0.5)
        pl.axhline(pi0,color='g',linewidth=2,label='$\pi$')
        pl.axhline(1,color='r',linestyle='--',linewidth=2,label='null')
        pl.legend()
        pl.xlabel('permutation p-value')
        pl.ylabel('density')
        pl.ylim((0,10))
        pl.savefig(ouFprefix + '.' + p_ID + '.pdf')
        pl.close('all')

fdrPlot('QTL-lm-permutation-Cmp')


def fdrCmp(ouFprefix):
    for ip, p_ID in enumerate(dataset.phenotype_ID[0:P_max:2]):
        lm = qtl.test_lm(snps=geno[sample_idx],pheno=phenotype_vals.values)
        qv, qv_info = fdr.qvalues(lm.pvalues[1])
        pl.figure(figsize=[8,8])
        pl.plot((lm.pvalues[1]),(qv),'.',label = 'qv')
        pl.plot((lm.pvalues[1]),(lm.pvalues[1]*geno.shape[1]),'.',label = 'Bonferroni')
        pl.plot(lm.pvalues[1],lm.pvalues[1],'k--')
        pl.yscale("log")
        pl.xscale("log")
        pl.xlabel('p-value')
        pl.ylabel('adjuste p-value')
        pl.legend(loc='lower right')
        pl.savefig(ouFprefix + '.' + p_ID + '.pdf')
        pl.close('all')
fdrCmp('QTL-lm-fdr-Cmp')

