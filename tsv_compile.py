import glob
import sys
import math


def average(l):
    return sum(l) / float(len(l))

quant = {}
for fyle in glob.glob('*.tsv'):
    counts,tpm,i = [],[],0
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            counts.append(float(line.strip().split()[3]))
            tpm.append(float(line.strip().split()[4]))
            i+=1
    quant[fyle.strip().split(".")[0]] = [counts,tpm]


print i

with open(sys.argv[1]+'_counts.csv','w') as outp:
    outp.write('WT_1,WT_2,WT_3,dA_1,dA_2,dA_3,dG_1,dG_2,dG_3,dAG_1,dAG_2,dAG_3\n')
    q = 0
    while q <= i-1:
        print
        outp.write(str(quant['WT_1'][0][q])+",")
        outp.write(str(quant['WT_2'][0][q])+",")
        outp.write(str(quant['WT_3'][0][q])+",")
        outp.write(str(quant['dA_1'][0][q])+",")
        outp.write(str(quant['dA_2'][0][q])+",")
        outp.write(str(quant['dA_3'][0][q])+",")
        outp.write(str(quant['dG_1'][0][q])+",")
        outp.write(str(quant['dG_2'][0][q])+",")
        outp.write(str(quant['dG_3'][0][q])+",")
        outp.write(str(quant['dAdG_1'][0][q])+",")
        outp.write(str(quant['dAdG_2'][0][q])+",")
        outp.write(str(quant['dAdG_3'][0][q])+"\n")
        q+=1

with open(sys.argv[1]+'_tpm.csv','w') as outp:
    outp.write('WT_1,WT_2,WT_3,dA_1,dA_2,dA_3,dG_1,dG_2,dG_3,dAG_1,dAG_2,dAG_3\n')
    q = 0
    while q <= i-1:
        outp.write(str(quant['WT_1'][1][q])+",")
        outp.write(str(quant['WT_2'][1][q])+",")
        outp.write(str(quant['WT_3'][1][q])+",")
        outp.write(str(quant['dA_1'][1][q])+",")
        outp.write(str(quant['dA_2'][1][q])+",")
        outp.write(str(quant['dA_3'][1][q])+",")
        outp.write(str(quant['dG_1'][1][q])+",")
        outp.write(str(quant['dG_2'][1][q])+",")
        outp.write(str(quant['dG_3'][1][q])+",")
        outp.write(str(quant['dAdG_1'][1][q])+",")
        outp.write(str(quant['dAdG_2'][1][q])+",")
        outp.write(str(quant['dAdG_3'][1][q])+"\n")
        q+=1
