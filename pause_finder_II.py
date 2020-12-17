from __future__ import division
import sys
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import collections
import numpy as np
import matplotlib.pyplot as plt

def cruncher(lyst):
    return round((lyst.count('T')+lyst.count('A'))/len(lyst)*100,3)

def seq_strip(fyle):
    seqs,crunched = [],[]
    for record in SeqIO.parse(fyle, "fasta"):
        seqs.append(str(record.seq))
    q = 0
    while q < len(seqs[0])-1:
        crunched.append(cruncher([item[q] for item in seqs]))
        q +=1
    return crunched

final = {}
for fyle in glob.glob('*.fa'):
    final[fyle] = seq_strip(fyle)

with open(sys.argv[1],'w') as outp:
    outp.write('Pos,Type,Value\n')
    q = 0
    while q <= len(final['G.fa'])-1:
        outp.write(str(q+1)+",G,"+str(final['G.fa'][q])+"\n")
        outp.write(str(q+1)+",AG,"+str(final['AG.fa'][q])+"\n")
        outp.write(str(q+1)+",A,"+str(final['A.fa'][q])+"\n")
        q+=1
