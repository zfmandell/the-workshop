#!/usr/bin/env python

import glob
import sys
import argparse

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i

def collect_infomation():
    #Snags a directory worth of CT files,runs terminal report on each CT file, outputs list of terminal reports
    information = []
    for fyle in glob.glob('*.ct'):
        information.append(terminal_report(fyle))
    return information

def terminal_report(fyle):
    #generates comprehensive information list of lists, indexed as mentioned above (#1-5)
    return_list = []
    len_fold = file_len(fyle)
    with open(fyle,'r') as inp:
        trash = inp.readline()

    return return_list


def main():
    parser = argparse.ArgumentParser(description='will create .csv file of all delta values at all called peaks ')
    parser.add_argument('outfile',type=str,default=None,help=' base name of output files')
    args = parser.parse_args()


    with open(args.outfile,'w') as outp:

if __name__ == '__main__':
    main()
