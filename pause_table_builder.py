#!/usr/bin/env python

from __future__ import division
import sys
import itertools
import math

# order == WT,dA,dG,dAG #

WT = sys.argv[2]
mutant = sys.argv[3]

pops_WT = {}
with open(sys.argv[2],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        if float(line.strip().split('\t')[8]) < 0:
            pops_WT[int(line.strip().split('\t')[1])] = ['-',abs(float(line.strip().split('\t')[8])),abs(float(line.strip().split('\t')[9]))]
        else:
            pops_WT[int(line.strip().split('\t')[1])] = ['+',abs(float(line.strip().split('\t')[8])),abs(float(line.strip().split('\t')[9]))]

pops_mutant = {}
with open(sys.argv[3],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        if float(line.strip().split('\t')[8]) < 0:
            pops_mutant[int(line.strip().split('\t')[1])] = ['-',abs(float(line.strip().split('\t')[8])),abs(float(line.strip().split('\t')[9]))]
        else:
            pops_mutant[int(line.strip().split('\t')[1])] = ['+',abs(float(line.strip().split('\t')[8])),abs(float(line.strip().split('\t')[9]))]

final = []
for key,value in pops_WT.iteritems():
    if key in pops_mutant.keys():
        final.append([key,pops_WT[key][0],pops_WT[key][1],pops_WT[key][2],(pops_WT[key][1]/pops_mutant[key][1]),(pops_WT[key][2]/pops_mutant[key][2])])
    else:
        final.append([key,pops_WT[key][0],pops_WT[key][1],pops_WT[key][2],'NA','NA'])

with open(sys.argv[1],'w') as outp:
    outp.write('coord,strand,WT_Counts,WT_Score,FC_Counts,FC_Score\n')
    outp.write('\n'.join(','.join(map(str, x)) for x in final))
