def g012tohdf5(g012_file):
    genotype = hdf.create_group('genotype')
    col_header = genotype.create_group('col_header')
    row_header = genotype.create_group('row_header')
    #load position and meta information
    indv_file = g012_file + '.indv'
    pos_file  = g012_file + '.pos'
    sample_ID = sp.loadtxt(indv_file,dtype='str')
    pos  = sp.loadtxt(pos_file,dtype='str')
    chrom = pos[:,0]
    pos   = sp.array(pos[:,1],dtype='int')

    row_header.create_dataset(name='sample_ID',data=sample_ID)
    col_header.create_dataset(name='chrom',data=chrom)
    col_header.create_dataset(name='pos',data=pos)
    M = sp.loadtxt(g012_file,dtype='uint8')
    snps = M[:,1::]
    genotype.create_dataset(name='matrix',data=snps,chunks=(snps.shape[0],min(10000,snps.shape[1])),compression='gzip')

g012tohdf5('Yeast-Genotype-SA.012')
