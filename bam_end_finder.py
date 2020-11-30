#!/usr/bin/env python

from __future__ import division
import argparse
import collections

def median(lst):
    n = len(lst)
    s = sorted(lst)
    return (sum(s[n//2-1:n//2+1])/2.0, s[n//2])[n % 2] if n else None

def peak_gen(fyle):
    output = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            output[int(line.strip().split(',')[0])] = line.strip().split(',')[1]
    return output

def end_gen(fyle):
    genome,output = {},[]
    for item in xrange(1,5000000):
        genome[item] = 0
    with open(fyle,'r') as inp:
        for line in inp:
            if len(line.strip().split('\t')) == 15:
                genome[int(line.strip().split('\t')[3])] += 1
    for k, v in collections.OrderedDict(sorted(genome.items())).iteritems():
        output.append(v)
    return output

def peak_placement(peaks,rev):
    output = {}
    for key,value in peaks.iteritems():
        if value == '-':
            portion = rev[key-2:key+3]
            try:
                output[key] = list(reversed([x / median(portion) for x in portion]))
            except ZeroDivisionError:
                pass
        else:
                pass
    return output

def outpwrite(placed,fyle):
    with open(fyle,'w') as outp:
        outp.write('pos,value\n')
        for key,value in placed.iteritems():
            outp.write('-2,'+str(value[0])+'\n')
            outp.write('-1,'+str(value[1])+'\n')
            outp.write('0,'+str(value[2])+'\n')
            outp.write('1,'+str(value[3])+'\n')
            outp.write('2,'+str(value[4])+'\n')

def main():
    parser = argparse.ArgumentParser(description='finds the distribution of ends around a terminator')
    parser.add_argument('-fwd',type=str,help='fwd strand cov file')
    parser.add_argument('-rev',type=str,help='rev strand cov file')
    parser.add_argument('-peaks',type=str,help='flat text file of terminators')
    parser.add_argument('-outfile',type=str,help='name of output file')
    args = parser.parse_args()

    peaks = peak_gen(args.peaks)

    #fwd_cov = end_gen(args.fwd)
    rev_cov = end_gen(args.rev)

    placed = peak_placement(peaks,rev_cov)

    outpwrite(placed,args.outfile)

if __name__ == '__main__':
    main()
