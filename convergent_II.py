from __future__ import division
import sys

gaps = []
with open(sys.argv[1],'r') as inp:
    for line in inp:
        gaps.append([int(line.strip().split()[0]),int(line.strip().split()[1])])

q,i,ends = 0,0,[]
with open(sys.argv[2],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        for item in gaps:
            if item[0] < int(line.strip()) < item[1]:
                q += 1
                ends.append(int(line.strip()))
        i+=1
with open(sys.argv[2],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        if int(line.strip()) not in ends:
            print line.strip()

print i
print q
print round((q/i)*100,3)
