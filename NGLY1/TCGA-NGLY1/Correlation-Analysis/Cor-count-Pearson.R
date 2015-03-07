load('../UCEC.rda')
library(DESeq2)
#rld=rlog(dds)
rld=dds
rld.assay=assay(rld)
rld.cor=cor(t(rld.assay), use="complete.obs")
ngly1=rld.cor[grepl('NGLY1',rownames(rld.cor)),]
ngly1.plus.0.5=ngly1[which(ngly1>=(0.5))]
ngly1.plus.0.6=ngly1[which(ngly1>=(0.6))]
ngly1.minus.0.6=ngly1[which(ngly1<=(-0.6))]
ngly1.minus.0.5=ngly1[which(ngly1<=(-0.5))]
write.table(as.data.frame(ngly1.minus.0.5),file='NGLY1-count-Pearson-minus-0.5',quote=F,col.names=F)
write.table(as.data.frame(ngly1.minus.0.6),file='NGLY1-count-Pearson-minus-0.6',quote=F,col.names=F)
write.table(as.data.frame(ngly1.plus.0.6),file='NGLY1-count-Pearson-plus-0.6',quote=F,col.names=F)
write.table(as.data.frame(ngly1.plus.0.5),file='NGLY1-count-Pearson-plus-0.5',quote=F,col.names=F)
