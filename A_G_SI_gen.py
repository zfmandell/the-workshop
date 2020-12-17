import sys

A,G,AG,SI = [],[],[],[]
with open(sys.argv[1],'r') as inp:
    firstline = inp.readline()
    for line in inp:

        gene = line.strip().split(',')[2].split('+')[0]

        A_term = float(line.strip().split(',')[5])
        G_term = float(line.strip().split(',')[7])
        WT = float(line.strip().split(',')[3])

        if A_term >= 25 and G_term < 10:
            A.append(str(gene))
        elif G_term >= 25 and A_term < 10:
            G.append(str(gene))
        elif A_term < 10 and G_term < 10:
            if WT > 70:
                SI.append(str(gene))
        elif A_term >= 25 and G_term >= 25:
            AG.append(str(gene))

with open('SI_genes.txt','w') as outp:
    for item in SI:
        outp.write(item+'\n')
with open('A_genes.txt','w') as outp:
    for item in A:
        outp.write(item+'\n')
with open('G_genes.txt','w') as outp:
    for item in G:
        outp.write(item+'\n')
with open('AG_genes.txt','w') as outp:
    for item in AG:
        outp.write(item+'\n')
