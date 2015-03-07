rt=read.table("nationwidechildrens.org_clinical_patient_ucec.txt", sep="\t", header=T)
ngly1=rt[rt[,1] %in% c("TCGA-D1-A17Q","TCGA-B5-A0JY","TCGA-D1-A103","TCGA-B5-A11N"),]
