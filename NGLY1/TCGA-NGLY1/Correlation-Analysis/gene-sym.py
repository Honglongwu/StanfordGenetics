import sys
inFile = open(sys.argv[1])
ouFile = open(sys.argv[1]+'.gene', 'w')
for line in inFile:
    fields = line.split('|')
    if fields[0] != '?':
        ouFile.write(fields[0] + '\n')
inFile.close()
ouFile.close()
