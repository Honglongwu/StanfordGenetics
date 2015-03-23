inFile = open('TCGA_RNA-Seq_GeneCounts-UCEC')
head = inFile.readline().split('\t')
for line in inFile:
    line = line.strip()
    fields = line.split('\t')
    gene = fields[0]
    count = [int(x) for x in fields[1:]]
    count_sorted = sorted(count)
    maxCount = count_sorted[-1]
    minCount = count_sorted[0]
    maxCount_i = count.index(maxCount)
    minCount_i = count.index(minCount)
    maxCount_s = head[maxCount_i+1]
    minCount_s = head[minCount_i+1]
    maxCount_second = count_sorted[-2]
    minCount_second = count_sorted[1]
    print('%s\t%s\t%s\t%s\t%s\t%s\t%s'%(gene,maxCount_s,minCount_s,maxCount,minCount,maxCount_second,minCount_second))
inFile.close()
