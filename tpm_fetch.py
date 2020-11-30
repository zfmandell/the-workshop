from __future__ import division
import glob
import sys
import math

def average(l):
    return sum(l) / float(len(l))


quant = {}
for fyle in glob.glob('*.tsv'):
    tpm = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            tpm[line.strip().split()[0]] = float(line.strip().split()[3])
    quant[fyle.strip().split(".")[0]] = tpm

mutant = sys.argv[1]

genes_WT,genes_mutant = {},{}
for key,value in quant.iteritems():
    if 'WT' == key.split('_')[0]:
        for key,sub_value in value.iteritems():
            if key in genes_WT.keys():
                genes_WT[key].append(sub_value)
            else:
                genes_WT[key] = [0,sub_value]
    if mutant == key.split('_')[0]:
        for key,sub_value in value.iteritems():
            if key in genes_mutant.keys():
                genes_mutant[key].append(sub_value)
            else:
                genes_mutant[key] = [0,sub_value]

final_WT = {}
for key,value in genes_WT.iteritems():
    final_WT[key] = average(value[1:])

final_mutant = {}
for key,value in genes_mutant.iteritems():
    final_mutant[key] = average(value[1:])

final = {}
for key,value in final_WT.iteritems():
    try:
        final[key] = value/final_mutant[key]
    except ZeroDivisionError:
        final[key] = 'NA'


with open(sys.argv[2],'w') as outp:
    outp.write('transcript,FC\n')
    for key,value in final.iteritems():
        if value != 0:
            try:
                outp.write(str(key)+','+str(math.log(value,2))+'\n')
            except TypeError:
                outp.write(str(key)+',NA\n')
        else:
            outp.write(str(key)+',NA\n')
