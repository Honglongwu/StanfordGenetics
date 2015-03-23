geneCounts = read.table('TCGA_RNA-Seq_GeneCounts',sep='\t',header=T,row.names=1)
sampleAnnot = read.table('TCGA_RNA-Seq_SampleAnnot',sep='\t',header=T,row.names=1)
save(geneCounts,sampleAnnot,file = 'TCGA-RNASeq.rda')
