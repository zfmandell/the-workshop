import sys
import glob

for fyle in glob.glob("./ends/*"):
    ends,total,sample = [],{},fyle.split(".")[0]
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            ends.append(line.strip())

    with open(sys.argv[1],'r') as inp:
        firstline = inp.readline()
        for line in inp:
            total[line.strip().split(",")[0]] = ",".join(line.split(",")[1:])

    with open(fyle+'_compiled.csv','w') as outp:
        outp.write(firstline)
        for item in sorted(ends):
            outp.write(str(item)+","+total[item])
