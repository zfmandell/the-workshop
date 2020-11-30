import sys

DA,DG,DADG,SI = [],[],[],[]
with open(sys.argv[1],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        A = float(line.strip().split(',')[1])
        G = float(line.strip().split(',')[2])
        if A >= 25 and -10 < G < 10:
            DA.append(line)
        elif G >= 25 and -10 < A < 10:
            DG.append(line)
        elif A >= 25 and G >= 25:
            DADG.append(line)
        elif -10 < G < 10 and -10 < A < 10:
            SI.append(line)
with open(sys.argv[2],'w') as outp:
    outp.write(firstline)
    for item in DA:
        outp.write(item)
