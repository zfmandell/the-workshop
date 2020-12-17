#!/usr/bin/env python

from __future__ import division
import argparse
from glob import glob

def main():
    parser = argparse.ArgumentParser(description='creates merged index file')
    parser.add_argument('outfile',type=str,default=None,help='name of output files')
    parser.add_argument('peaks',type=str,default=None,help='name of peak file')
    args = parser.parse_args()

    with open(args.peaks,'r') as inp:
        firstline = inp.readline()
        conditions,peaks = [str(firstline.strip().split(",")[0]),str(firstline.strip().split(",")[1])],[]
        secondline = inp.readline()
        for line in inp:
            peaks.append([int(line.strip().split(',')[0]),line.strip().split(',')[1]])

    fwd = [x for x in sorted(glob("./*fwd*")) if x.split("_")[1] in conditions]
    rev = [x for x in sorted(glob("./*rev*")) if x.split("_")[1] in conditions]

    with open(conditions[0]+"_"+conditions[1]+args.outfile,'w') as outp:
        outp.write('peak,')
        outp.write(','.join(["_".join([x.split("_")[1],x.split("_")[2]]) for x in fwd]))
        outp.write("\n")
        for item in peaks:
            to_add = []
            if item[1] == 'fwd':
                for fyle in fwd:
                    with open(fyle,'r') as fp:
                        for i, line in enumerate(fp):
                            if i == item[0]+1:
                                to_add.append(str(line.strip().split("\t")[3]))
                            elif i > item[0]+1:
                                break

                outp.write(str(item[0])+",")
                outp.write(",".join(to_add))
                outp.write("\n")

            elif item[1] == 'rev':
                for fyle in rev:
                    with open(fyle,'r') as fp:
                        for i, line in enumerate(fp):
                            if i == item[0]+1:
                                to_add.append(str(line.strip().split("\t")[3]))
                            elif i > item[0]+1:
                                break

                outp.write(str(item[0])+",")
                outp.write(",".join(to_add))
                outp.write("\n")

if __name__ == '__main__':
    main()
