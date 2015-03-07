inFile = open('TCGA-ESP-count-enriched-total.filtered.fdr')
ouFile = open('TCGA-ESP-count-enriched-total.filtered.fdr.formated', 'w')
ouFile.write('\t'.join(['gene','cancerType','mutationType','No.Mutation(Cancer)','No.Cancer','No.Mutation(Normal)','No.Normal','p-value','q-value']) + '\n')
for line in inFile:
    line = line.strip()
    fields = line.split('\t')
    gene = fields[0]
    if gene.find('_')==-1:
        ouFile.write('\t'.join([fields[0], fields[2], fields[6], fields[4], fields[5], fields[10], fields[11], fields[14], fields[15]]) + '\n')
inFile.close()
ouFile.close()
