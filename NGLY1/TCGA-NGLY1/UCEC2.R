library(DESeq2)
load('UCEC.rda')
### DESeq2
#ngly1.control.tmp = sapply(ngly1.control, as.integer)
#rownames(ngly1.control.tmp)=rownames(ngly1.control)
#ngly1.control=ngly1.control.tmp
dds=DESeqDataSetFromMatrix(countData=ngly1.control, colData=ngly1.control.annotation,design=~sample+condition)
dds <- DESeq(dds)
dds.results=results(dds)
dds.results=dds.results[order(dds.results$padj),]
dds.results.significant=dds.results[which(dds.results$padj<0.05),]
write.table(dds.results.significant, 'UCEC-3NGLY1-36Control-Significant-Genes-Factor-Sample',quote=F,col.names=NA)

#rld = rlog(dds, blind=FALSE)
#
#pdf("plot_PCA-NGLY1-Control.pdf")
#print(plotPCA(rld))
#dev.off()
#
#pdf("plot_MA-NGLY1-Control.pdf")
#plotMA(dds.results, alpha=0.01)
#dev.off()
