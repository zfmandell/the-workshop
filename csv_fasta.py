import sys
from Bio import SeqIO

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases

csv = sys.argv[1]
fasta = sys.argv[2]
upstream = int(sys.argv[3])
downstream = int(sys.argv[4])
output = sys.argv[5]

record = SeqIO.read(fasta, "fasta")
reference = str(record.seq)

coords = {}
with open(csv,'r') as inp:
    firstline = inp.readline()
    for line in inp:
        coords[int(line.strip().split(',')[0])] = line.strip().split(',')[1]

with open(output,'w') as outp:
    for key,value in coords.items():
        if value == '+':
            outp.write('>'+str(key)+'_'+value+'\n')
            outp.write(reference[key-upstream:key+downstream]+'\n')
        else:
            outp.write('>'+str(key)+'_'+value+'\n')
            sequence = reference[key-downstream:key+upstream]
            outp.write(reverse_complement(sequence)+'\n')
