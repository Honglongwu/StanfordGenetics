rt = read.table('UCEC-NGLY1-ga.txt', sep='\t')
rt[grepl('TCGA-D1-A17Q-01', unname(unlist(rt[1,])))]
rt[grepl('TCGA-B5-A0JY-01', unname(unlist(rt[1,])))]
rt[grepl('TCGA-D1-A103-01', unname(unlist(rt[1,])))]
rt[grepl('TCGA-B5-A11N-01', unname(unlist(rt[1,])))]

select=seq(3,dim(rt)[2],2)
count=unname(unlist(rt[2,select]))
count= as.numeric(as.character(count))
sample=unname(unlist(rt[1,select]))
df = data.frame(sample, count)
save(count,sample,file='UCEC-NGLY1-ga.rda')


library(ggplot2)
pdf('UCEC-NGLY1-ga-count-distribution-bin50.pdf')
color = rep(c("white", "red", "white", "red", "white"), times=c(5,3,12,1,53))
ggplot(df, aes(x=count)) + geom_histogram(binwidth=50,colour="black", fill=color) + ylab("Number of Samples") + xlab('NGLY1 Expression (RNA-Seq reads count)')
#ggplot(df, aes(x=count)) + geom_histogram(colour="black", fill="white") + ylab("Number of Samples") + xlab('NGLY1 Expression (RNA-Seq raw reads count)')
dev.off()

pdf('UCEC-NGLY1-ga-count-distribution-bin100.pdf')
color = rep(c("white", "red", "white", "red", "white"), times=c(3,2,6,1,27))
ggplot(df, aes(x=count)) + geom_histogram(binwidth=100,colour="black",fill=color) + ylab("Number of Samples") + xlab('NGLY1 Expression (RNA-Seq reads count)')
dev.off()

pdf('UCEC-NGLY1-ga-count-distribution-bin100-2.pdf')
color = rep(c("white", "red", "white","blue", "white"), times=c(3,2,1,1,32))
ggplot(df, aes(x=count)) + geom_histogram(binwidth=100,colour="black",fill=color) + ylab("Number of Samples") + xlab('NGLY1 Expression (RNA-Seq reads count)')
dev.off()

