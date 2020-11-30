#!/usr/bin/env python

from __future__ import division
from glob import glob
import sys
from Bio import SeqIO
import re

def chunks(seq, window_size):
    it = iter(seq)
    win = [it.next() for cnt in xrange(window_size)] # First window
    yield win
    for e in it: # Subsequent windows
        win[:-1] = win[1:]
        win[-1] = e
        yield win

def reverse_complement(seq):
    alt_map = {'ins':'0'}
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    for k,v in alt_map.iteritems():
        seq = seq.replace(k,v)
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.iteritems():
        bases = bases.replace(v,k)
    return bases

terms = {}
with open(sys.argv[1],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        terms[int(line.strip().split(',')[0])] = line.strip().split(',')[1]

for record in SeqIO.parse(sys.argv[2], "fasta"):
    genome = str(record.seq)

competitors = {}
for fyle in glob('*.ct'):
    POT = int(fyle.strip().split('_')[0])
    for item in terms.keys():
        if item-2 <= POT <= item+2:
            with open(fyle,'r') as inp:
                temp_paired = []
                temp_nucleotide = []
                firstline = inp.readline()
                for line in inp:
                    temp_paired.append(int(line.strip().split()[4]))
                    temp_nucleotide.append(str(line.strip().split()[1]))
                    temp_chunks = chunks(temp_paired,3)
                q = 0
                for sub_item in temp_chunks:
                    if sub_item == [0,0,0]:
                        try:
                            break
                        except ZeroDivisionError:
                            break
                    else:
                        q+=1
                strand = terms[item]
                stem = ''.join(temp_nucleotide[:q])
                if strand == '+':
                    upstream = genome[POT-250:POT].replace('T','U')
                    location = upstream.rfind(stem)
                    if location != -1:
                        competition = upstream[location-50:location+int(round(len(stem)/2))].replace('T','U')
                        competitors[item] = competition
                else:
                    upstream = reverse_complement(genome[POT:POT+250]).replace('T','U')
                    location = upstream.rfind(stem)
                    if location != -1:
                        competition = upstream[location-50:location+int(round(len(stem)/2))].replace('T','U')
                        competitors[item] = competition

with open(sys.argv[3],'w') as outp:
    for key,value in competitors.iteritems():
        outp.write('>'+str(key)+'\n'+str(value)+'\n')
