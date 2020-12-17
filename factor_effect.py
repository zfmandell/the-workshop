import sys

dependent,strong_independent,weak_independent = [],[],[]
with open(sys.argv[1],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        if float(line.strip().split(",")[5]) >= 25.0:
            dependent.append(int(line.strip().split(",")[0]))
        elif float(line.strip().split(",")[5]) < 5.0 and float(line.strip().split(",")[3]) >= 70.0:
            strong_independent.append(int(line.strip().split(",")[0]))
        elif float(line.strip().split(",")[5]) < 5.0 and float(line.strip().split(",")[3]) <= 30:
            weak_independent.append(int(line.strip().split(",")[0]))

with open('D_'+sys.argv[1].split("_")[1],'w') as outp:
    with open(sys.argv[1],'r') as inp:
        firstline = inp.readline()
        outp.write(firstline)
        for line in inp:
            if int(line.strip().split(",")[0]) in dependent:
                outp.write(line)

with open('SI_'+sys.argv[1].split("_")[1],'w') as outp:
    with open(sys.argv[1],'r') as inp:
        firstline = inp.readline()
        outp.write(firstline)
        for line in inp:
            if int(line.strip().split(",")[0]) in strong_independent:
                outp.write(line)

with open('WI_'+sys.argv[1].split("_")[1],'w') as outp:
    with open(sys.argv[1],'r') as inp:
        firstline = inp.readline()
        outp.write(firstline)
        for line in inp:
            if int(line.strip().split(",")[0]) in weak_independent:
                outp.write(line)
