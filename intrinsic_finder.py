from __future__ import division
import sys
from scipy import stats

terms = {}
with open(sys.argv[1],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        terms[int(line.strip().split(',')[0])] = [line.strip().split(',')[1],line.strip().split(',')[2]]

final = {}
with open(sys.argv[2],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        end = int(line.strip().split(',')[0])
        strand = line.strip().split(',')[1]

        for key,value in terms.iteritems():
            if key-3 <= end <= key+3 and value[0] == strand:
                final[key] = delta


with open(sys.argv[3],'w') as outp:
    outp.write('POT,strand\n')
    for key,value in final.iteritems():
        outp.write(str(key)+','+value+'\n')"""
