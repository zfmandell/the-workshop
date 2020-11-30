
import sys
from glob import glob
from Bio import SeqIO

reference = {}
with open(sys.argv[1],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        coord = line.strip().split(',')[0]
        strand = line.strip().split(',')[1]
        cv = line.strip().split(',')[2]
        pT = line.strip().split(',')[5]
        rP = line.strip().split(',')[6]
        reference[str(coord)+'_'+str(strand)] = [str(cv),str(pT),str(rP)]

fastas = {}
with open(sys.argv[2], "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        coord = str(record.id).split('strand')[0].split(':')[1]
        strand = str(record.id).split(':')[2].split('_')[0]
        fastas[str(coord)+'_'+str(strand)] = str(record.seq)

deltas = {}
for fyle in glob(str(sys.argv[3])+'/*.ct'):
    print(fyle)
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        if len(firstline.split()) == 5:
            coord = str(firstline.split()[4]).split('strand')[0].split(':')[1]
            strand = str(firstline.split()[4]).split(':')[2].split('_')[0]
            deltas[str(coord)+'_'+str(strand)] = str(firstline.split()[3])

for key in reference.keys():
    if key in deltas.keys():
        to_append = [fastas[key],deltas[key]]
        for item in to_append:
            reference[key].append(item)
    else:
        to_append = [fastas[key],str(0)]
        for item in to_append:
            reference[key].append(item)

with open(sys.argv[4],'w') as outp:
    outp.write('end,Strand,Cv,upstream 50 nt,dG,%T,Relative Position\n')
    for key,value in reference.items():
        coord = key.split('_')[0]
        strand = key.split('_')[1]
        cv = value[0]
        nuc = value[3]
        dG = value[4]
        pT = value[1]
        rP = value[2]
        outp.write(coord+','+strand+','+cv+','+nuc+','+dG+','+pT+','+rP+'\n')
