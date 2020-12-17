import sys

ends = []
with open(sys.argv[1],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        ends.append(int(line.strip()))

with open(sys.argv[2],'r') as inp:
    firstline = inp.readline()
    with open(sys.argv[3],'w') as outp:
        outp.write(firstline)
        for line in inp:
            if int(line.strip().split(",")[0]) in ends:
                outp.write(line)
