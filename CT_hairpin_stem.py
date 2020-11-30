#!/usr/bin/env python

from __future__ import division
from glob import glob
import sys
from collections import Counter

def chunks(seq, window_size):
    it = iter(seq)
    win = [it.next() for cnt in xrange(window_size)] # First window
    yield win
    for e in it: # Subsequent windows
        win[:-1] = win[1:]
        win[-1] = e
        yield win

final_GC = {}
final_length_loop = {}


terms = []
with open(sys.argv[1],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        terms.append(int(line.strip().split(',')[0]))

for fyle in glob('*.ct'):
    POT = int(fyle.strip().split('_')[0])
    for item in terms:
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
                            coun = Counter(temp_nucleotide[0:q])
                            final_GC[item] = ((coun['G']+coun['C'])/len(temp_nucleotide[0:q]))*100
                            break
                        except ZeroDivisionError:
                            break
                    else:
                        q+=1


                final_length_loop[item] = ''.join(temp_nucleotide[:q])



"""
with open(sys.argv[2],'w') as outp:
    outp.write('POT,GC\n')
    for key,value in final_GC.iteritems():
        outp.write(str(key)+','+str(value)+'\n')
"""

with open(sys.argv[2],'w') as outp:
    outp.write('POT,percA,percC,percG,percU\n')
    for key,value in final_length_loop.iteritems():
        outp.write('>'+str(key)+'\n'+str(value)+'\n')
