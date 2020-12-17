#!/usr/bin/env python


from __future__ import division
import sys
from glob import glob

def read_wig(wig_fyle,lyst):
    #create list of coverage values, position dependent
    return_dict = {}
    with open(wig_fyle,"r") as inp:
        for line in inp:
            if int(line.strip().split('\t')[0]) in lyst:
                try:
                    return_dict[int(line.strip().split('\t')[0])] = abs(int(line.strip().split('\t')[1])/float(line.strip().split('\t')[3]))
                except ZeroDivisionError:
                    return_dict[int(line.strip().split('\t')[0])] = 'NA'
    return return_dict

wild,mutant = [],[]
with open(sys.argv[1],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        wt = float(line.strip().split(',')[8])
        mut = float(line.strip().split(',')[9])

        if wt < 1 :
            wild.append(int(line.strip().split(',')[0]))
        elif mut < 1:
            mutant.append(int(line.strip().split(',')[0]))

print mutant

wt_final,mutant_final = {},{}
w,m = {},{}

for fyle in glob('*.wig'):
    wt_final[fyle] = read_wig(fyle,wild)
    mutant_final[fyle] = read_wig(fyle,mutant)


for key,value in wt_final.iteritems():
    for subkey,subvalue in value.iteritems():
        if subkey in w.keys():
            w[subkey].append(subvalue)
        else:
            w[subkey] = [subvalue]

for key,value in mutant_final.iteritems():
    for subkey,subvalue in value.iteritems():
        if subkey in m.keys():
            w[subkey].append(subvalue)
        else:
            w[subkey] = [subvalue]

print wt_final

with open(sys.argv[2],'w') as outp:
    for key,value in w.iteritems():
        outp.write(str(key)+',')
        outp.write(','.join(map(str,value))+'\n')

with open(sys.argv[3],'w') as outp:
    for key,value in m.iteritems():
        outp.write(str(key)+',')
        outp.write(','.join(map(str,value))+'\n')

#print len(wild), ' = # holes in WT'
#print len(mutant), ' = # holes in Muant'
