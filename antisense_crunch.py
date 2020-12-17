#!/usr/bin/env python

from __future__ import division
import argparse
from glob import glob
from scipy import stats

def read_gff3(gff3):
    #reads in fasta, returns dict of sequences, dict of strand info for gene
    return_dict = {}
    with open(gff3,'r') as inp:
        for line in inp:
            if line[0] != '#' and len(line.strip().split()) != 0:
                if str(line.strip().split()[2]) == 'gene':
                    return_dict[(line.strip().split()[8].split(";")[2].split("=")[1]).replace("-",".")] = [line.strip().split()[6],int(line.strip().split()[3]),int(line.strip().split()[4])]
    return return_dict


def read_coverage(wig_fyle):
    #create list of coverage values, position dependent
    return_list = []
    with open(wig_fyle,"r") as inp:
        for line in inp:
            return_list.append(int(line.split("\t")[2]))
    return return_list

def main():
    parser = argparse.ArgumentParser(description='creates merged index file')
    parser.add_argument('outfile',type=str,default=None,help='name of output files')
    parser.add_argument('gff',type=str,default=None,help='name of peak file')
    parser.add_argument('condition',type=str,default=None,help='condition of interest')
    args = parser.parse_args()

    fwd = [x for x in sorted(glob("./*fwd*")) if x.split("_")[0][2:] == args.condition]
    rev = [x for x in sorted(glob("./*rev*")) if x.split("_")[0][2:] == args.condition]
    total,peaks = {},read_gff3(args.gff)

    for fyle in fwd:
        cov_dict[str(fyle)] = read_coverage(fyle)
    for fyle in rev:
        cov_dict[str(fyle)] = read_coverage(fyle)

    for key,value in peaks.iteritems():
        if value[0] == '+':
            temp = []
            for fyle in rev:
                temp.append(sum(cov_dict[fyle][value[0]:value[1]+1]))
        else:
            temp = []
            for fyle in fwd:
                temp.append(sum(cov_dict[fyle][value[0]:value[1]+1]))
        total[key] = temp
    with open(args.outfile,'w') as outp:
        outp.write('unit,'+args.condition+'_1,'+args.condition+'_2,'+args.condition+'_3\n')
        for key,value in total.iteritems():
            outp.write(str(key)+","str(value[0])+","str(value[1])+","str(value[2])+"\n")

if __name__ == '__main__':
    main()
