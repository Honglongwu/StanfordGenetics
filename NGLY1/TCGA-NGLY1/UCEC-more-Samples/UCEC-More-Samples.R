###
load('UCEC.rda')
selection = NGLY1.count.t<=600 & NGLY1.count.t >= 500 & annotation$platform=='ga' & annotation$race=='WHITE' & annotation$age <= 70 & annotation$age >= 50
control=count[,selection]
control.annotation = annotation[selection,]
write.table(control,file="UCEC-36control",quote=F,col.names=NA)
###ngly1=count[,grepl("TCGA-D1-A17Q-01|TCGA-B5-A0JY-01|TCGA-D1-A103-01|TCGA-B5-A11N-01",colnames(count))]
ngly1=count[,grepl("TCGA-D1-A17Q-01|TCGA-B5-A0JY-01|TCGA-B5-A11N-01",colnames(count))]
ngly1.annotation = annotation[annotation$sample %in% colnames(ngly1),]
write.table(ngly1,file="UCEC-3NGLY1",quote=F,col.names=NA)

### merge ngly1 and control
condition = rep(c('ngly1','control'), times=c(dim(ngly1.annotation)[1], dim(control.annotation)[1]))
ngly1.control = merge(ngly1, control, by.x=0, by.y=0)
rownames(ngly1.control) = ngly1.control[,1]
ngly1.control=ngly1.control[2:dim(ngly1.control)[2]]
ngly1.control.annotation = cbind(rbind(ngly1.annotation, control.annotation),condition)
#ngly1.control = merge(ngly1, control, by.x=1, by.y=1)

### DESeq2
#ngly1.control.tmp = sapply(ngly1.control, as.integer)
#rownames(ngly1.control.tmp)=rownames(ngly1.control)
#ngly1.control=ngly1.control.tmp
library(DESeq2)
dds=DESeqDataSetFromMatrix(countData=ngly1.control, colData=ngly1.control.annotation,design=~condition)
dds <- DESeq(dds)
dds.results=results(dds)
dds.results=dds.results[order(dds.results$padj),]
dds.results.significant=dds.results[which(dds.results$padj<0.05),]
save(dds,ngly1.control,ngly1.control.annotation,file="UCEC-More.rda")

#rld = rlog(dds, blind=FALSE)
#
#pdf("plot_PCA-NGLY1-Control.pdf")
#print(plotPCA(rld))
#dev.off()
#
#pdf("plot_MA-NGLY1-Control.pdf")
#plotMA(dds.results, alpha=0.01)
#dev.off()
