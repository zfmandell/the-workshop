#!/usr/bin/env python

import argparse
import itertools
import operator
from glob import glob

def peak_clean(lyst):
    return_list = []
    set_list = sorted(set(lyst))
    for item in set_list:
        try:
            if item < set_list[set_list.index(item)+1]-3:
                return_list.append(item)
        except IndexError:
            pass
    return return_list

def main():
    parser = argparse.ArgumentParser(description='creates merged index file')
    parser.add_argument('outfile',type=str,default=None,help='name of output files')
    parser.add_argument('condition_A',type=str,default=None,help='biological condition A')
    parser.add_argument('condition_B',type=str,default=None,help='biological condition B')
    args = parser.parse_args()

    peaks_fwd,peaks_rev,conditions = [],[],[args.condition_A,args.condition_B]
    for fyle in glob("./*fwd*"):
        if fyle.split("_")[3] in conditions:
            with open(fyle,'r') as inp:
                for line in inp:
                    peaks_fwd.append(int(line.split("\t")[2]))

    for fyle in glob("./*rev*"):
        if fyle.split("_")[3] in conditions:
            with open(fyle,'r') as inp:
                for line in inp:
                    peaks_rev.append(int(line.split("\t")[2]))

    clean_fwd,clean_rev = peak_clean(peaks_fwd),peak_clean(peaks_rev)
    fwd_stranded,rev_stranded = [[x,'fwd'] for x in clean_fwd],[[x,'rev'] for x in clean_rev]
    all_peaks = sorted(list(itertools.chain.from_iterable([fwd_stranded,rev_stranded])),key=operator.itemgetter(0))

    with open(args.outfile,'w') as outp:
        outp.write('peak,strand\n')
        for item in all_peaks:
            outp.write(str(item[0])+","+str(item[1])+"\n")

if __name__ == '__main__':
    main()
