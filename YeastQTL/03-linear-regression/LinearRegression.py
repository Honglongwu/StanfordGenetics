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

geno_reader = gr.genotype_reader_tables('Yeast-Genotype.hdf5')
pheno_reader = phr.pheno_reader_tables('Yeast-Phenotype.hdf5')
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

'''
for ip, p_ID in enumerate(phenotype_ID):
    pl.figure()
    plot_normal(phenotype_vals.values[:,ip],alpha=0.8)
    pl.title("%s" % p_ID)
    pl.savefig('Phenotype-%s.pdf'%p_ID)

phenotype_vals_boxcox, maxlog = preprocess.boxcox(phenotype_vals.values)
phenotype_vals_ranks = preprocess.rankStandardizeNormal(phenotype_vals.values)

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
'''

lm = qtl.test_lm(snps=geno[sample_idx],pheno=phenotype_vals.values)
pvalues_lm = pd.DataFrame(data=lm.pvalues.T,index=dataset.geno_ID,columns=phenotype_ID)
ouFile = open('QTL-lm-significant','w')
for i in range(pvalues_lm.shape[1]):
    n = -1
    for x in pvalues_lm.ix[:,i]:
        n += 1
        if x < 0.000001:
            ouFile.write("%s\t%s_%s\t%s\t%s\t%s"%(n,pos['chrom'][n], pos['pos'][n],i*2,pvalues_lm.columns[i],x) + '\n')
ouFile.close()

'''
for ip, p_ID in enumerate(phenotype_ID):
    pl.figure(figsize=[12,4])
    plot_manhattan(posCum=pos['pos_cum'],pv=pvalues_lm[p_ID].values,chromBounds=chromBounds,thr_plotting=0.05)
    pl.title(p_ID)
    pl.savefig('QTL-untransformed-manhatton-%s.pdf'%p_ID)

lm = qtl.test_lm(snps=geno[sample_idx],pheno=phenotype_vals_boxcox)
pvalues_lm = pd.DataFrame(data=lm.pvalues.T,index=dataset.geno_ID,columns=phenotype_ID)
for ip, p_ID in enumerate(phenotype_ID):
    pl.figure(figsize=[12,4])
    plot_manhattan(posCum=pos['pos_cum'],pv=pvalues_lm[p_ID].values,chromBounds=chromBounds,thr_plotting=0.05)
    pl.title(p_ID)
    pl.savefig('QTL-boxcox-manhatton-%s.pdf'%p_ID)

lm = qtl.test_lm(snps=geno[sample_idx],pheno=phenotype_vals_ranks)
pvalues_lm = pd.DataFrame(data=lm.pvalues.T,index=dataset.geno_ID,columns=phenotype_ID)
for ip, p_ID in enumerate(phenotype_ID):
    pl.figure(figsize=[12,4])
    plot_manhattan(posCum=pos['pos_cum'],pv=pvalues_lm[p_ID].values,chromBounds=chromBounds,thr_plotting=0.05)
    pl.title(p_ID)
    pl.savefig('QTL-ranks-manhatton-%s.pdf'%p_ID)
'''

phenotype_ID = dataset.phenotype_ID[72:73]
for ip, p_ID in enumerate(phenotype_ID):
    pl.figure(figsize=[12,4])
    plot_manhattan(posCum=pos['pos_cum'],pv=pvalues_lm[p_ID].values,chromBounds=chromBounds,thr_plotting=0.05)
    pl.title(p_ID)
    pl.savefig('QTL-untransformed-manhatton-%s.pdf'%p_ID)
 
 
    
        



