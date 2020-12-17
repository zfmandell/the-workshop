#!/usr/bin/env python

import argparse
import numpy as np
from scipy.stats import ttest_ind
import operator

def translator(lyst,base):
    for item in base:
        if item[0]-3 <= lyst[0] <= item[0]+3 and lyst[1] == item[1][0]:
            return item[0]
            break
    return None

def WT_term(fyle):
    return_dict = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            if line.strip().split(",")[1] == 'fwd':
                return_dict[int(line.strip().split(",")[0])] = ['+',line.strip().split(",")[2],line.strip().split(",")[3]]
            else:
                return_dict[int(line.strip().split(",")[0])] = ['-',line.strip().split(",")[2],line.strip().split(",")[3]]
    return return_dict

def WT_rep(fyle):
    return_dict = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            if line.strip().split(",")[1] == 'fwd':
                return_dict[int(line.strip().split(",")[0])] = ['+',float(line.strip().split(",")[2]),float(line.strip().split(",")[3]),float(line.strip().split(",")[4])]
            else:
                return_dict[int(line.strip().split(",")[0])] = ['-',float(line.strip().split(",")[2]),float(line.strip().split(",")[3]),float(line.strip().split(",")[4])]
    return return_dict

def Exp_term(fyle,WT):
    return_dict = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            if line.strip().split(",")[1] == 'fwd':
                if str(translator([int(line.strip().split(",")[0]),'+'],WT.items())) != 'None':
                    return_dict[translator([int(line.strip().split(",")[0]),'+'],WT.items())] = ['+',line.strip().split(",")[2],line.strip().split(",")[3]]
            else:
                if str(translator([int(line.strip().split(",")[0]),'-'],WT.items())) != 'None':
                    return_dict[translator([int(line.strip().split(",")[0]),'-'],WT.items())] = ['-',line.strip().split(",")[2],line.strip().split(",")[3]]
    return return_dict

def Exp_rep(fyle,WT):
    return_dict = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            if line.strip().split(",")[1] == 'fwd':
                if translator([int(line.strip().split(",")[0]),'+'],WT.items()) != None:
                    return_dict[translator([int(line.strip().split(",")[0]),'+'],WT.items())] = ['+',float(line.strip().split(",")[2]),float(line.strip().split(",")[3]),float(line.strip().split(",")[4])]
            else:
                if translator([int(line.strip().split(",")[0]),'-'],WT.items()) != None:
                    return_dict[translator([int(line.strip().split(",")[0]),'-'],WT.items())] = ['-',float(line.strip().split(",")[2]),float(line.strip().split(",")[3]),float(line.strip().split(",")[4])]
    return return_dict

def term_build(fyle):
    return_dict = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            return_dict[int(line.strip().split(",")[0])] = line.strip().split(",")[1:]
    return return_dict

def mean(l):
    return reduce(lambda x, y: x + y, l) / len(l)

def stats(a,b):
    if str(ttest_ind(a,b)[1]) == 'nan':
        return 'ND'
    else:
        return format_e(ttest_ind(a,b)[1])

def format_e(n):
    a = '%E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('WildType',type=str,help='')
    parser.add_argument('Experimental',type=str,help='')
    parser.add_argument('intrinsic_key',type=str,help='')
    parser.add_argument('outfile',type=str,help='')
    args = parser.parse_args()

    WT_terms,WT_reps = WT_term(args.WildType+'_term_eff.csv'),WT_rep(args.WildType+'_term_eff_rep.csv')
    Exp_terms,Exp_reps = Exp_term(args.Experimental+'_term_eff.csv',WT_terms),Exp_rep(args.Experimental+'_term_eff_rep.csv',WT_terms)
    intrinsic = term_build(args.intrinsic_key)

    with open(args.outfile,'w') as outp:
        outp.write('POT,Strand,Relative Position,T ('+args.WildType+'),T ('+args.Experimental+'),T,p-value,upstream sequence,hairpin sequence,U-tract,G\n')
        for item in sorted(intrinsic.items(), key=operator.itemgetter(0)):
            if translator([item[0],item[1][0]],WT_terms.items()) in WT_terms.keys():
                outp.write(str(item[0])+',')
                outp.write(str(item[1][0])+',')
                outp.write(str(item[1][1])+',')
                outp.write(WT_terms[translator([item[0],item[1][0]],WT_terms.items())][1]+",")
                if translator([item[0],item[1][0]],WT_terms.items()) in Exp_terms.keys():
                    outp.write(Exp_terms[translator([item[0],item[1][0]],WT_terms.items())][1]+",")
                    outp.write(str(float(WT_terms[translator([item[0],item[1][0]],WT_terms.items())][1])-float(Exp_terms[translator([item[0],item[1][0]],WT_terms.items())][1]))+",")
                    outp.write(str(stats(WT_reps[translator([item[0],item[1][0]],WT_terms.items())][1:],Exp_reps[translator([item[0],item[1][0]],WT_terms.items())][1:])+","))
                else:
                    outp.write('0,')
                    outp.write(WT_terms[translator([item[0],item[1][0]],WT_terms.items())][1]+",")
                    outp.write(stats(WT_reps[translator([item[0],item[1][0]],WT_terms.items())][1:],[0.0,0.0,0.0])+",")
                outp.write(str(item[1][2])+',')
                outp.write(str(item[1][3])+',')
                outp.write(str(item[1][4])+',')
                outp.write(str(item[1][5])+'\n')

if __name__ == '__main__':
    main()
