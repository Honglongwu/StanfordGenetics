import matplotlib
matplotlib.use('Agg')
import limix.io.genotype_reader as gr
import limix.io.phenotype_reader as phr
import scipy as sp
import pylab as pl
import limix.stats.pca as pca
from limix.stats.geno_summary import *

def genotypeQC(inF):
    prefix = inF.split('.hdf5')[0]

    genotype_file = gr.genotype_reader_tables(inF)
    genotype = genotype_file.getGenotypes()
    pos=genotype_file.getPos()
    
    K = sp.dot(genotype,genotype.T)
    K2 = sp.dot(genotype[:,0:2],genotype[:,0:2].T)
    pl.figure()
    pl.subplot(121)
    pl.imshow(K,aspect='auto')
    pl.subplot(122)
    pl.imshow(K2,aspect='auto')
    pl.savefig(prefix + "-dotplot.pdf")
    
    pl.figure()
    af = calc_AF(genotype,minor=1)
    pl.hist(af['af'],50,normed=True)
    pl.ylabel('density')
    pl.title('Yeast')
    pl.xlabel('allele frequency')
    pl.savefig(prefix + "-af.pdf")
    
    
    pl.figure()
    maf = 10.**(-sp.arange(6))
    N = [(af['af']>=m).sum() for m in maf]
    pl.plot(maf,N,'k-')
    pl.xscale('log')
    pl.xlabel('maf cutoff')
    pl.ylabel('number of variants')
    pl.title('Yeast')
    pl.savefig(prefix + "-maf.pdf")
    
    
    pl.figure()
    ld = calc_LD(genotype,pos['pos'].values)
    pl.plot(ld[0],ld[1])
    pl.ylabel('r^2')
    pl.xlabel('dist (bp)')
    pl.title('Yeast')
    pl.savefig(prefix + "-ld.pdf")
    
    pl.figure()
    [snpsX,snpsW] = pca.PCA(genotype[:,0:5000],components=10)
    pl.plot(snpsX[:,0],snpsX[:,1],'k.')
    pl.xlabel('PC1')
    pl.ylabel('PC2')
    pl.title('Yeast')
    pl.savefig(prefix + "-pca.pdf")

genotypeQC('Yeast-Genotype-noMissing-WA.hdf5')
