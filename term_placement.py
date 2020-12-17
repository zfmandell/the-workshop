#!/usr/bin/env python

import argparse
import operator

def build_end_dict(fyle):
    return_dyct = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            if line.strip().split(",")[1] == 'fwd':
                return_dyct[int(line.strip().split(",")[0])] = '+'
            else:
                return_dyct[int(line.strip().split(",")[0])] = '-'
    return return_dyct

def read_gff3(gff3):
    #reads in fasta, returns dict of sequences, dict of strand info for gene
    fwd,rev = {},{}
    with open(gff3,'r') as inp:
        for line in inp:
            if line[0] != '#' and len(line.strip().split()) != 0:
                if str(line.strip().split()[2]) == 'gene' and str(line.strip().split()[6]) == '+':
                    fwd[(line.strip().split()[8].split(";")[2].split("=")[1]).replace("-",".")] = [int(line.strip().split()[3]),int(line.strip().split()[4])]
                elif str(line.strip().split()[2]) == 'gene' and str(line.strip().split()[6]) == '-':
                    rev[(line.strip().split()[8].split(";")[2].split("=")[1]).replace("-",".")] = [int(line.strip().split()[3]),int(line.strip().split()[4])]
    return sorted(fwd.items(), key=operator.itemgetter(1)),sorted(rev.items(), key=operator.itemgetter(1))

def end_place_fwd(end,gff):
    for item in gff:
        index = gff.index(item)
        if end > item[1][1] and end < gff[gff.index(item)+1][1][0]:
            behind = end-item[1][1]
            fwd = gff[gff.index(item)+1][1][0]-end
            if behind < fwd:
                return str(item[0])+'+'+str(behind)+'nt'
            else:
                return str(gff[gff.index(item)+1][0])+'-'+str(fwd)+'nt'

def end_place_rev(end,gff):
    for item in gff:
        index = gff.index(item)
        if end > item[1][1] and end < gff[gff.index(item)+1][1][0]:
            fwd = end-item[1][1]
            behind = gff[gff.index(item)+1][1][0]-end
            if behind < fwd:
                return str(item[0])+'+'+str(behind)+'nt'
            else:
                return str(gff[gff.index(item)+1][0])+'-'+str(fwd)+'nt'

def misc_place(ends,place_dict):
    misc = [i for i in ends.keys() if place_dict[i] == None]

    for item in misc:
         place_dict[item] = 'misc'

    return place_dict

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('gff3',type=str,default=None,help='')
    parser.add_argument('ends',type=str,default=None,help='')
    parser.add_argument('outfile',type=str,default=None,help='name of output file')
    args = parser.parse_args()

    fwd,rev = read_gff3(args.gff3)
    ends = build_end_dict(args.ends)

    final = {}
    for key,value in ends.iteritems():
        if value == '+':
            final[key] = end_place_fwd(key,fwd)
        else:
            final[key] = end_place_rev(key,rev)

    final = misc_place(ends,final)

    with open(args.outfile,'w') as outp:
        outp.write('POT,strand,annotation\n')
        for item in sorted(final.items(), key=operator.itemgetter(0)):
            outp.write(str(item[0])+","+ends[item[0]]+","+item[1]+"\n")

if __name__ == '__main__':
    main()
