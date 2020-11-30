#!/usr/bin/env python

import glob
import os
import argparse
import subprocess

def end_dict(dyct,fyle):
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            if line.strip().split(",")[1] == 'fwd':
                dyct[int(line.strip().split(",")[0])] = '+'
            else:
                dyct[int(line.strip().split(",")[0])] = '-'

def term_list(lyst,fyle):
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            lyst.append([int(line.strip().split(",")[0]),line.strip().split(",")[1]])

def read_gff3(gff3):
    #reads in fasta, returns dict of sequences, dict of strand info for gene
    return_list = []
    with open(gff3,'r') as inp:
        for line in inp:
            if line[0] != '#' and len(line.strip().split()) != 0:
                if str(line.strip().split()[2]) == 'gene':
                    return_list.append([int(line.strip().split()[3]),int(line.strip().split()[4]),line.strip().split()[6]])
    return return_list

def gene_check(end,genes,dyct):
    for item in genes:
        if (item[0] < end < item[1]) and dyct[end] == item[2]:
            return True
    return False

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('all_ends',type=str,help='list of ends')
    parser.add_argument('intrinsic_confirmed',type=str,help='list of confirmed intrinsic ends')
    parser.add_argument('gff3',type=str,help='genbank for BSub')
    args = parser.parse_args()

    ends,intrinsic,genes = {},[],read_gff3(args.gff3)
    term_list(intrinsic,args.intrinsic_confirmed)
    end_dict(ends,args.all_ends)

    os.chdir('PS')
    names = glob.glob('*.ps')
    os.mkdir("../final_ps",0755)
    for item in names:
        if int(item.split("_")[0].split(":")[1]) not in [x[0] for x in intrinsic] and gene_check(int(item.split("_")[0].split(":")[1]),genes,ends) == False:
            subprocess.call('mv '+str(item)+' ../final_ps',shell=True)

if __name__ == '__main__':
    main()
