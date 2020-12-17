#!/usr/bin/env python

import sys
import operator

genes = []
with open(sys.argv[1],'r') as inp:
    for line in inp:
        if len(line.strip().split('\t')) > 1:
            if line.strip().split('\t')[2] == 'CDS':
                gene = line.strip().split('\t')[8].split('=')[7].split(';')[0]
                start = int(line.strip().split('\t')[3])
                end = int(line.strip().split('\t')[4])
                strand = line.strip().split('\t')[6]
                genes.append([gene,start,end,strand])


annotated_genes = {}
ends = {}
with open(sys.argv[2],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        ends[int(line.strip().split(',')[0])] = line.strip().split(',')[1]

for end,strand in ends.iteritems():
    for item in genes:
        if end <= item[2] and end >= item[1]:
            if strand == item[3]:
                annotated_genes[end] = item[0]+'_CDS'

genes_plus = [x for x in genes if x[3] == '+']
genes_minus = [x for x in genes if x[3] == '-']

for end,strand in ends.iteritems():
    print end
    if end in annotated_genes.keys():
        pass
    else:
        if strand == '+':
            for item in genes_plus:
                try:
                    if end > item[2] and end < genes_plus[genes_plus.index(item)+1][1]:
                        annotated_genes[end] = item[0]+'+'+str(end-item[2])
                except IndexError:
                    annotated_genes[end] = 'misc.'

        else:
            try:
                for item in genes_minus:
                    if end > item[2] and end < genes_minus[genes_minus.index(item)+1][1]:
                        annotated_genes[end] = genes_minus[genes_minus.index(item)+1][0]+'+'+str(genes_minus[genes_minus.index(item)+1][1]-end)
            except IndexError:
                annotated_genes[end] = 'misc.'

for item in ends.keys():
    if item not in annotated_genes.keys():
        annotated_genes[item] = 'misc.'

for item in ends.keys():
    if item not in annotated_genes.keys():
        print item

with open(sys.argv[3],'w') as outp:
    outp.write('POS,Relative Position\n')
    for key,value in annotated_genes.iteritems():
        outp.write(str(key)+','+str(value)+'\n')
