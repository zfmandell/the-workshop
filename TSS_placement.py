import sys
from Bio import SeqIO

query,split,TSS = {},{},[]

for record in SeqIO.parse(sys.argv[1], "fasta"):
    query['NCBI'] = str(record.seq)

with open(sys.argv[2],'r') as inp:
    for line in inp:
        TSS.append(int(line.strip().split('\t')[3]))

for item in TSS:
    if item < 50 or len(query['NCBI'])-item < 50:
        pass
    else:
        split[str(item)] = query['NCBI'][item-50:item+51]

with open(sys.argv[3],'w') as outp:
    for key,value in split.iteritems():
        outp.write('>'+key+'\n'+value+'\n')
