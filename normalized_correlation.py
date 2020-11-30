#!/usr/bin/env python

import argparse
from collections import Counter
from itertools import islice
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from scipy.stats.stats import pearsonr

def read_in_rtsc(rtsc_fyle):
    '''Reads a <.rtsc> file into a dictionary'''
    information = {}
    with open(rtsc_fyle,'r') as f:
        while True:
            next_n_lines = list(islice(f, 3))
            if not next_n_lines:
                break
            transcript,stops,empty_line = [n.strip() for n in next_n_lines]
            information[transcript] = [int(x) for x in stops.split('\t')]
    return information


def get_covered_transcripts(coverage_fyle):
    '''Reads in a standard overlapped coverage file'''
    info = {}
    with open(coverage_fyle,'r') as f:
        for line in f:
            info[line.strip()] = None
    return info

def filter_dictonary(in_dict,filter_dict):
    '''Removes entries from a dictionary'''
    for k in in_dict.keys():
        if k not in filter_dict:
            del in_dict[k]

def get_sum(alist):
    '''returns the sum of a list, removing NAs'''
    return sum([x for x in alist if x != 'NA'])

def normalize_rtsc(rtsc_dict):
    information = {}
    for key,value in rtsc_dict.iteritems():
        summed = get_sum(value)
        to_return = []
        for item in value:
            to_return.append(float(item)/summed)
        information[key] = to_return

    return information

def main():
    parser = argparse.ArgumentParser(description='Determines normalized correlation between two RTSC')
    parser.add_argument("rtsc_one",type=str,help="<.rtsc> file")
    parser.add_argument("rtsc_two",type=str,help="<.rtsc> file")
    parser.add_argument('restrict',default = None, help = '<.txt > Limit analysis to these specific transcripts')
    parser.add_argument('outname',type=str,help='name of the outfile')
    args = parser.parse_args()

    one,two = read_in_rtsc(args.rtsc_one),read_in_rtsc(args.rtsc_two)

    if args.restrict != None:
        covered = get_covered_transcripts(args.restrict)
        filter_dictonary(one,covered)
        filter_dictonary(two,covered)

    n_one,n_two = normalize_rtsc(one),normalize_rtsc(two)

    correl = {}
    for key, value in n_one.iteritems():
        correl[key] = (pearsonr(value,n_two[key]))


    with open(args.outname,"w") as outp:
        outp.write("transcript,correlation")
        outp.write("\n")
        for key, value in correl.iteritems():
            outp.write(str(key)+","+str(value[0])+"\n")

if __name__ == '__main__':
    main()
