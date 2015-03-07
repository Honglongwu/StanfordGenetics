import limix.io.genotype_reader as gr
import limix.io.phenotype_reader as phr
import scipy as sp
import pylab as pl
import limix.stats.pca as pca
from limix.stats.geno_summary import *


genotype_file = gr.genotype_reader_tables("Yeast-Genotype.hdf5")
genotype = genotype_file.getGenotypes()
pos=genotype_file.getPos()

K = sp.dot(genotype,genotype.T)
K2 = sp.dot(genotype[:,0:2],genotype[:,0:2].T)
pl.figure()
pl.subplot(121)
pl.imshow(K,aspect='auto')
pl.subplot(122)
pl.imshow(K2,aspect='auto')
pl.savefig("genotype-qc-dotplot.pdf")

pl.figure()
af = calc_AF(genotype,minor=2)
pl.hist(af['af'],50,normed=True)
pl.ylabel('density')
pl.title('Yeast')
pl.xlabel('allele frequency')
pl.savefig("genotype-qc-af.pdf")


pl.figure()
maf = 10.**(-sp.arange(6))
N = [(af['af']>=m).sum() for m in maf]
pl.plot(maf,N,'k-')
pl.xscale('log')
pl.xlabel('maf cutoff')
pl.ylabel('number of variants')
pl.title('Yeast')
pl.savefig("genotype-qc-maf.pdf")


pl.figure()
ld = calc_LD(genotype,pos['pos'].values)
pl.plot(ld[0],ld[1])
pl.ylabel('r^2')
pl.xlabel('dist (bp)')
pl.title('Yeast')
pl.savefig("genotype-qc-ld.pdf")

pl.figure()
[snpsX,snpsW] = pca.PCA(genotype[:,0:5000],components=10)
pl.plot(snpsX[:,0],snpsX[:,1],'k.')
pl.xlabel('PC1')
pl.ylabel('PC2')
pl.title('Yeast')
pl.savefig("genotype-qc-pca.pdf")
