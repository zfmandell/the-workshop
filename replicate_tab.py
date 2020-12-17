#!/usr/bin/env python

from __future__ import division
import sys
from scipy.stats import spearmanr
import numpy as np

def read_quant(quant_fyle):
    #create list of coverage values, position dependent
    return_list = []
    with open(quant_fyle,"r") as inp:
        firstline = inp.readline()
        for line in inp:
            return_list.append(float(line.strip().split('\t')[3]))
    return np.array(return_list)

first_wig,second_wig,third_wig = read_quant(sys.argv[1]),read_quant(sys.argv[2]),read_quant(sys.argv[3])

one_two_final = spearmanr(first_wig,second_wig)[0]
one_three_final = spearmanr(first_wig,third_wig)[0]
two_three_final = spearmanr(second_wig,third_wig)[0]

with open(sys.argv[4],'w') as outp:
    outp.write('file 1,file 2,# conserved,spearman r\n')
    outp.write(str(sys.argv[1])+','+str(sys.argv[2])+','+str(one_two_final)+'\n')
    outp.write(str(sys.argv[1])+','+str(sys.argv[3])+','+str(one_three_final)+'\n')
    outp.write(str(sys.argv[2])+','+str(sys.argv[3])+','+str(two_three_final))
