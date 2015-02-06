data = read.table("TCGA-ESP-count-enriched-total.filtered", sep="\t")
#fdr=p.adjust(data[,15],"fdr")
bh=p.adjust(data[,15],"BH")
#data$fdr = fdr
data$bh = bh
data = data[order(data$bh),]
write.table(data, file = "TCGA-ESP-count-enriched-total.filtered.fdr", quote=F, col.names=F,row.names=F,sep="\t")
