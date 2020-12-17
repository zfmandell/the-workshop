#!/usr/bin/env python

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import random

fasta_sequences,fasta_dict =SeqIO.parse(open(sys.argv[1]),'fasta'),{}
for fasta in fasta_sequences:
    fasta_dict[str(fasta.id)] = str(fasta.seq)

final = []
with open(sys.argv[2],'w') as outp:
    while len(final) < int(sys.argv[3]):
        rando = random.choice(list(fasta_dict.keys()))
        if rando not in final:
            final.append(rando)

    for item in final:
        outp.write('>'+str(item)+"\n"+str(fasta_dict[item])+'\n')
