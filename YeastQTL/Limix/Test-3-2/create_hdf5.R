library(rhdf5)

X=data.frame(s1=c('AA','AB','BC'),s2=c('AC','BC','BB'),s3=c('BC','BB','AB'))
Y=data.frame(p1=c(1,2,2),s2=c(3,2,2),s3=c(1,1,2))
myfile='test.hdf5'

  # X is the genotype matrix, which has individuals in columns and genomic positions in rows.
  # Y is the phenotype matrix, which has individuals in rows and different phenotypes in columns.
  
  # Note:
  # 1. X must have sample_IDs as colnames and genomic positions (in the format "chr1-9999") as rownames
  # 2. Y must have pheno_IDs as colnames
  # 3. The columns of X and rows of Y must be in the same order (i.e. correspond to same indiduals)
  splitted.info = strsplit(rownames(X), "-")
  npos = nrow(X)
  alleles = matrix(NA, 2, npos)
  chrom = as.numeric(substr(sapply(splitted.info, function(x)x[1]), 4, 5))
  pos = as.numeric(sapply(splitted.info, function(x)x[2]))
  genotype = X
  sample_ID = colnames(X)
  phenotype_ID = colnames(Y)
  phenotype = as.matrix(Y)
  colnames(genotype) <- rownames(genotype) <- NULL
  colnames(phenotype) <- rownames(phenotype) <- NULL
  
  
  h5createFile(myfile)
  # Genotype
  h5createGroup(myfile, "genotype")
  h5createGroup(myfile, "genotype/col_header")
  h5write(alleles, file=myfile, name="genotype/col_header/alleles")
  h5write(chrom, file=myfile, name="genotype/col_header/chrom")
  h5write(pos, file=myfile, name="genotype/col_header/pos")
  h5write(genotype, file=myfile, name="genotype/matrix")
  h5createGroup(myfile, "genotype/row_header")
  h5write(sample_ID, file=myfile, name="genotype/row_header/sample_ID")
  # Phenotype
  h5createGroup(myfile, "phenotype")
  h5createGroup(myfile, "phenotype/col_header")
  h5write(phenotype_ID, file=myfile, name="phenotype/col_header/phenotype_ID")
  h5write(t(phenotype), file=myfile, name="phenotype/matrix")
  h5createGroup(myfile, "phenotype/row_header")
  h5write(sample_ID, file=myfile, name="phenotype/row_header/sample_ID")
  H5close()
