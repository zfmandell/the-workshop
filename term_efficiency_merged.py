#!/usr/bin/env python

from __future__ import division
import argparse
import glob
from scipy.stats import ttest_ind
import numpy as np

def median(lst):
    n = len(lst)
    s = sorted(lst)
    return (sum(s[n//2-1:n//2+1])/2.0, s[n//2])[n % 2] if n else None

def find_plus(coverage_list,peak):
    if peak < 11:
        window = 11 - peak
        return coverage_list[:peak]+coverage_list[len(coverage_list)-window-1:]
    else:
        return coverage_list[peak-11:peak]

def find_minus(coverage_list,peak):
    if peak > len(coverage_list)-11:
        window = len(coverage_list)-peak
        return coverage_list[peak+1:]+coverage_list[:12-window]
    else:
        return coverage_list[peak+1:peak+12]

def read_coverage(cov_fyle):
    #create list of coverage values, position dependent
    return_list = []
    with open(cov_fyle,"r") as inp:
        for line in inp:
            return_list.append(int(float((line.strip().split('\t')[2]))))
    return return_list

def term_eff(up,down):
    try:
        eff = (median(up)-median(down))/median(up)*100
        if eff > 0:
            return eff
        else:
            return 0
    except ZeroDivisionError:
        return 0


def main():
    parser = argparse.ArgumentParser(description='creates all termination efficiency files')
    parser.add_argument('-peaks',type=str,default=None,help='name of peak file')
    args = parser.parse_args()

    peaks,hashable_terms = [],[]
    with open(args.peaks,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            hashable_terms.append(line.strip().split(',')[0]+'_'+line.strip().split(',')[1])

    WT,dA,dG,dAdG = [],[],[],[]
    WT_fwd,dA_fwd,dG_fwd,dAdG_fwd = [],[],[],[]
    WT_rev,dA_rev,dG_rev,dAdG_rev = [],[],[],[]
    strains = [WT,dA,dG,dAdG]
    for fyle in glob.glob("*.cov"):
        if fyle.split('_')[0] == 'WT':
            WT.append(fyle)
            if fyle.split('_')[3] == 'fwd.cov':
                if len(WT_fwd) > 1:
                    temp_cov = read_coverage(fyle)
                    temp_list = [a + b for a, b in zip(WT_fwd,temp_cov)]
                    WT_fwd = temp_list
                else:
                    WT_fwd = read_coverage(fyle)
            if fyle.split('_')[3] == 'rev.cov':
                if len(WT_rev) > 1:
                    temp_cov = read_coverage(fyle)
                    temp_list = [a + b for a, b in zip(WT_rev,temp_cov)]
                    WT_rev = temp_list
                else:
                    WT_rev = read_coverage(fyle)
        elif fyle.split('_')[0] == 'dA':
            dA.append(fyle)
            if fyle.split('_')[3] == 'fwd.cov':
                if len(dA_fwd) > 1:
                    temp_cov = read_coverage(fyle)
                    temp_list = [a + b for a, b in zip(dA_fwd,temp_cov)]
                    dA_fwd = temp_list
                else:
                    dA_fwd = read_coverage(fyle)
            if fyle.split('_')[3] == 'rev.cov':
                if len(dA_rev) > 1:
                    temp_cov = read_coverage(fyle)
                    temp_list = [a + b for a, b in zip(dA_rev,temp_cov)]
                    dA_rev = temp_list
                else:
                    dA_rev = read_coverage(fyle)
        elif fyle.split('_')[0] == 'dG':
            dG.append(fyle)
            if fyle.split('_')[3] == 'fwd.cov':
                if len(dG_fwd) > 1:
                    temp_cov = read_coverage(fyle)
                    temp_list = [a + b for a, b in zip(dG_fwd,temp_cov)]
                    dG_fwd = temp_list
                else:
                    dG_fwd = read_coverage(fyle)
            if fyle.split('_')[3] == 'rev.cov':
                if len(dG_rev) > 1:
                    temp_cov = read_coverage(fyle)
                    temp_list = [a + b for a, b in zip(dG_rev,temp_cov)]
                    dG_rev = temp_list
                else:
                    dG_rev = read_coverage(fyle)
        elif fyle.split('_')[0] == 'dAdG':
            dAdG.append(fyle)
            if fyle.split('_')[3] == 'fwd.cov':
                if len(dAdG_fwd) > 1:
                    temp_cov = read_coverage(fyle)
                    temp_list = [a + b for a, b in zip(dAdG_fwd,temp_cov)]
                    dAdG_fwd = temp_list
                else:
                    dAdG_fwd = read_coverage(fyle)
            if fyle.split('_')[3] == 'rev.cov':
                if len(dAdG_rev) > 1:
                    temp_cov = read_coverage(fyle)
                    temp_list = [a + b for a, b in zip(dAdG_rev,temp_cov)]
                    dAdG_rev = temp_list
                else:
                    dAdG_rev = read_coverage(fyle)

    WT_final,dA_final,dG_final,dAdG_final = {},{},{},{}
    for item in hashable_terms:
        pos,strand = int(item.split('_')[0]),item.split('_')[1]
        if strand == '+':
            WT_final[item] = term_eff(find_plus(WT_fwd,pos),find_minus(WT_fwd,pos))
            dA_final[item] = term_eff(find_plus(dA_fwd,pos),find_minus(dA_fwd,pos))
            dG_final[item] = term_eff(find_plus(dG_fwd,pos),find_minus(dG_fwd,pos))
            dAdG_final[item] = term_eff(find_plus(dAdG_fwd,pos),find_minus(dAdG_fwd,pos))
        else:
            WT_final[item] = term_eff(find_minus(WT_rev,pos),find_plus(WT_rev,pos))
            dA_final[item] = term_eff(find_minus(dA_rev,pos),find_plus(dA_rev,pos))
            dG_final[item] = term_eff(find_minus(dG_rev,pos),find_plus(dG_rev,pos))
            dAdG_final[item] = term_eff(find_minus(dAdG_rev,pos),find_plus(dAdG_rev,pos))

    WT_dA,WT_dG,WT_dAdG = {},{},{}
    for item in hashable_terms:
        WT_dA[item] = [item.split('_')[0],item.split('_')[1],WT_final[item],dA_final[item],WT_final[item]-dA_final[item]]
        WT_dG[item] = [item.split('_')[0],item.split('_')[1],WT_final[item],dG_final[item],WT_final[item]-dG_final[item]]
        WT_dAdG[item] = [item.split('_')[0],item.split('_')[1],WT_final[item],dAdG_final[item],WT_final[item]-dAdG_final[item]]

    with open('WT_dA_total.csv','w') as outp:
        outp.write('peak,strand,Term-WT,Term-dA,dTerm\n')
        for key,value in WT_dA.iteritems():
            outp.write(','.join(map(str,value)))
            outp.write('\n')

    with open('WT_dG_total.csv','w') as outp:
        outp.write('peak,strand,Term-WT,Term-dG,dTerm\n')
        for key,value in WT_dG.iteritems():
            outp.write(','.join(map(str,value)))
            outp.write('\n')

    with open('WT_dAdG_total.csv','w') as outp:
        outp.write('peak,strand,Term-WT,Term-dAdG,dTerm\n')
        for key,value in WT_dAdG.iteritems():
            outp.write(','.join(map(str,value)))
            outp.write('\n')

if __name__ == '__main__':
    main()
