library(DESeq2)

load('TCGA-RNASeq.rda')
condition=rep('control',dim(sampleAnnot)[1])
sampleAnnot$condition=condition
C=read.table('NGLY1-Mutated-Samples')[,1]
sampleAnnot[sampleAnnot$sample %in% C,]$condition='mutated'

cancer = sampleAnnot[sampleAnnot$sample %in% C,]$cancer
select1 = sampleAnnot$cancer %in% cancer
select2 = sampleAnnot$sample_type == 'tumor'

wh = select1 & select2

count = geneCounts[,wh]
annot = droplevels(sampleAnnot[wh,])
annot$condition = factor(annot$condition,levels=c('control','mutated'))

dds = DESeqDataSetFromMatrix(count, annot, design=~cancer+condition)
dds = DESeq(dds)

res = results(dds)
res = res[order(res$padj), ]
res.sig = res[which(res$padj<0.05 & res$baseMean > 20),]
write.table( res, file="deNGLY1_2.txt", quote = FALSE, sep = "\t",  row.names = T, col.names=NA)
write.table( res.sig, file="deNGLY1_sig_2.txt", quote = FALSE, sep = "\t",  row.names = T, col.names=NA)
