import sys

genes = []
with open(sys.argv[1],'r') as inp:
    for line in inp:
        genes.append(line.strip())

final = []
with open(sys.argv[2],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        if line.strip().split(',')[0].split('-')[0].split('.')[1] in genes:
            final.append(str(line.strip().split(',')[1]))

with open(sys.argv[3],'w') as outp:
    outp.write('TPM\n')
    for item in final:
        outp.write(item+'\n')
