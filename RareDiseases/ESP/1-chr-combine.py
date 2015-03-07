import os

'''
#### not sorted by chr1, chr2 ...
Fs = os.listdir('.')
Fs.sort()
ouFile = open('ESP-snp-indels', 'w')

for F in Fs:
    if F.find('ESP6500SI-V2-SSA137.GRCh38-liftover') != -1:
        print(F)
        inFile = open(F)
        for line in inFile:
            if line[0] != '#':
                ouFile.write(line)
        inFile.close()
ouFile.close()

'''

chs = ['chr' + str(x) for x in range(1,23)] + ['chrX', 'chrY'] 
ouFile = open('ESP-snp-indels', 'w')
for ch in chs:
    inFile = open('ESP6500SI-V2-SSA137.GRCh38-liftover.' + ch + '.snps_indels.txt')
    for line in inFile:
        if line[0] != '#':
            ouFile.write(line)
    inFile.close()
ouFile.close()
