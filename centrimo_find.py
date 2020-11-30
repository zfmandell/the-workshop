import sys

sites = []
with open(sys.argv[1],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        if int(line.split()[2]) >= 255:
            sites.append(int(line.split()[1].split("_")[0].split(":")[1]))

to_keep = []
with open(sys.argv[2],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        if int(line.split(",")[0]) in sites:
            to_keep.append(line)

with open(sys.argv[3],'w') as outp:
    outp.write(firstline)
    for item in to_keep:
        outp.write(item)
