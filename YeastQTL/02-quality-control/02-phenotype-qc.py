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

geno_reader = gr.genotype_reader_tables('Yeast-Genotype-noMissing-WA.hdf5')
pheno_reader = phr.pheno_reader_tables('Yeast-Phenotype.hdf5')
dataset = data.QTLData(geno_reader=geno_reader,pheno_reader=pheno_reader)
geno = dataset.getGenotypes()
position = dataset.getPos()
pos,chromBounds = data_util.estCumPos(position=position,offset=0)
#P_max = len(dataset.phenotype_var_ID)
P_max = 6
phenotype_var_ID = dataset.phenotype_ID[1:P_max:2]
phenotype_vals, sample_idx = dataset.getPhenotypes(phenotype_var_ID)
N = geno.shape[0]
S = geno.shape[1]
P = phenotype_vals.shape[1]

def phenoPlot(phenotype_var_ID, phenotype_vals):
    for ip, p_ID in enumerate(phenotype_var_ID):
        pl.figure()
        plot_normal(phenotype_vals.values[:,ip],alpha=0.8)
        pl.title("%s" % p_ID)
        pl.savefig('Phenotype-%s.pdf'%p_ID)

phenoPlot(phenotype_var_ID, phenotype_vals)


