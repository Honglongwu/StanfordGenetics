inFile = open('morbidmap')
ouFile = open('omim-congenital-hereditary-inherited', 'w')
ouFile1 = open('omim-congenital', 'w')
ouFile2 = open('omim-hereditary', 'w')
ouFile3 = open('omim-inherited', 'w')

for line in inFile:
    line = line.strip()
    line_lower = line.lower()
    if line_lower.find('congenital') != -1 or line_lower.find('inherit') != -1 or line_lower.find('hereditary') != -1: 
        ouFile.write(line + '\n')
    if line_lower.find('congenital') != -1: 
        ouFile1.write(line + '\n')
    if line_lower.find('hereditary') != -1: 
        ouFile2.write(line + '\n')
    if line_lower.find('inherit') != -1: 
        ouFile3.write(line + '\n')

inFile.close()
ouFile1.close()
ouFile2.close()
ouFile3.close()
