
import sys
import numpy as np
import collections


def valCount(lst):
    res = {}
    for v in lst:
        try:
            res[v] += 1
        except KeyError:
            res[v] = 1
    return res

bedgraph = []
with open(sys.argv[1],'r') as inp:
    for line in inp:
        if float(line.strip().split('\t')[3]) < 0:
            bedgraph.append(str(line.strip().split('\t')[2])+'_'+'-')
        else:
            bedgraph.append(str(line.strip().split('\t')[2])+'_'+'+')

final,total = [],[]
with open(sys.argv[2],'rU') as inp:
    firstline = inp.readline()
    for line in inp:
        coord = str(line.strip().split(',')[0])
        strand = str(line.strip().split(',')[1])
        code = coord+'_'+strand
        for item in bedgraph:
            subcoord = int(item.split('_')[0])
            substrand = item.split('_')[1]
            if (subcoord-3 <= int(coord) <= subcoord+3) and (strand == substrand):
                final.append(code)
        total.append(code)

#u = [ x for x,y in valCount(total).iteritems() if y > 1 ]


print len(final)
print len(total)
if len(set(total)) == len(total):
    print 'true'
else:
    print 'false'
print np.setdiff1d(total,final)
