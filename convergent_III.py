#!/usr/bin/env python

from __future__ import division
import argparse
from glob import glob
from scipy import stats
import numpy as np

def find_plus(coverage_list,peak):
    if peak < 51:
        window = 51 - peak
        return coverage_list[:peak]+coverage_list[len(coverage_list)-window-1:]
    else:
        return coverage_list[peak-51:peak]

def find_minus(coverage_list,peak):
    if peak > len(coverage_list)-51:
        window = len(coverage_list)-peak
        return coverage_list[peak+1:]+coverage_list[:52-window]
    else:
        return coverage_list[peak+1:peak+52]

def read_coverage(wig_fyle):
    #create list of coverage values, position dependent
    return_list = []
    with open(wig_fyle,"r") as inp:
        for line in inp:
            return_list.append(int(line.split("\t")[2]))
    return return_list

def median(lst):
    n = len(lst)
    if n < 1:
            return None
    if n % 2 == 1:
            return sorted(lst)[n//2]
    else:
            return sum(sorted(lst)[n//2-1:n//2+1])/2.0

def main():
    parser = argparse.ArgumentParser(description='creates merged index file')
    parser.add_argument('peaks',type=str,default=None,help='name of peak file')
    parser.add_argument('fwd',type=str,default=None,help='condition of interest')
    parser.add_argument('rev',type=str,default=None,help='condition of interest')
    parser.add_argument('output',type=str,default=None,help='')
    args = parser.parse_args()

    fwd = [x for x in sorted(glob("./fwd/*fwd*")) if x.split("_")[0][2:] == args.condition]
    rev = [x for x in sorted(glob("./rev/*rev*")) if x.split("_")[0][2:] == args.condition]
    peaks,terms_total,cov_dict,fux = [],[],{},0

    with open(args.peaks,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            peaks.append([int(line.strip().split(',')[0]),line.strip().split(',')[1]])


    cov_dict[args.fwd] = read_coverage(args.fwd)
    cov_dict[args.rev] = read_coverage(args.rev)

    total = []
    for item in peaks:
        if item[1] == '+':
            up,down = find_plus(cov_dict[args.fwd],item[0]),find_minus(cov_dict[args.rev],item[0])
            if median(up) == 0 and median(down) != 0:
                total.append(np.log2(1/median(down)))
            elif median(down) == 0 and median(up) != 0:
                total.append(np.log2(median(up)/1))
            elif median(down) == 0 and median(up) == 0:
                total.append(np.log2(1/1))
            else:
                total.append(np.log2(median(up)/median(down)))
        if item[1] == '-':
            up,down= find_plus(cov_dict[args.fwd],item[0]),find_minus(cov_dict[args.rev],item[0])
            if median(up) == 0 and median(down) != 0:
                total.append(np.log2(median(down)/1))
            elif median(down) == 0 and median(up) != 0:
                total.append(np.log2(1/median(up)))
            elif median(down) == 0 and median(up) == 0:
                total.append(np.log2(1/1))
            else:
                total.append(np.log2(median(down)/median(up)))
    with open(args.output,'w') as outp:
        for item in total:
            outp.write(str(item)+"\n")


if __name__ == '__main__':
    main()
