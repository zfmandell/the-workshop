import sys
import numpy as np

strong = []
with open(sys.argv[1],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        strong.append(int(line.strip()))

independent = []
with open(sys.argv[2],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        independent.append(int(line.strip()))

with open('strong_'+sys.argv[2],'w') as outp:
    shared = [item for item in independent if item in strong]
    outp.write(firstline)
    for item in shared:
        outp.write(str(item)+"\n")

with open('weak_'+sys.argv[2],'w') as outp:
    not_shared = [item for item in independent if item not in strong]
    outp.write(firstline)
    for item in not_shared:
        outp.write(str(item)+"\n")
