#!/usr/bin/env python
import argparse

def term_build(lyst,fyle):
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            lyst.append([int(line.strip().split(",")[0]),line.strip().split(",")[1]])

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('my_term',type=str,help='my high term efficiency ends')
    parser.add_argument('jeet_term',type=str,help='Jeets high term efficiency ends')
    parser.add_argument('outfile',type=str,help='name of outfile')
    args = parser.parse_args()

    my_term,jeets_term,intrinsic = [],[],[]

    term_build(my_term,args.my_term)
    term_build(jeets_term,args.jeet_term)

    for item in my_term:
        if item[1] == 'fwd':
            item[1] = '+'
        else:
            item[1] = '-'

    for item in my_term:
        for sub_item in jeets_term:
            if (sub_item[0]-3 <= item[0] <= sub_item[0]+3) and item[1] == sub_item[1]:
                intrinsic.append(item)

    with open(args.outfile,'w') as outp:
        outp.write('ends,strand\n')
        for item in intrinsic:
            outp.write(str(item[0])+','+str(item[1])+"\n")

if __name__ == '__main__':
    main()
