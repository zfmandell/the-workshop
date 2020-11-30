#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def build_strand_list(fyle):
    return_list = []
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            return_list.append([int(line.strip().split(",")[0]),line.strip().split(",")[1]])
    return return_list

def build_strand_dict(fyle):
    return_dict = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            return_dict[int(line.strip().split(",")[0])] = line.strip().split(",")[1]
    return return_dict

def build_end_list(fyle):
    return_list = []
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            return_list.append(int(line.strip('\n')))
    return return_list

def build_duplex_dict(fyle):
    return_dict = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            return_dict[int(line.strip().split(",")[0])] = [line.strip().split(",")[1],line.strip().split(",")[2]]
    return return_dict

def duplex_translator(dyct,base):
    return_dict = {}
    for key,value in dyct.iteritems():
        for sub_item in base:
            if sub_item[0]-3 <= key <= sub_item[0]+3 and value[0] == sub_item[1]:
                return_dict[sub_item[0]] = value[1]
    return return_dict

def translator(lyst,base):
    for item in base:
        if item[0]-3 <= lyst[0] <= item[0]+3 and lyst[1] == item[1]:
            return item[0]

def deltaT_dict_build(fyle):
    return_dict = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            return_dict[int(line.strip().split(",")[0])] = float(line.strip().split(",")[2])
    return return_dict

def supplemental_build(fyle,base):
    return_dict = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            return_dict[translator([int(line.strip().split(",")[0]),line.strip().split(",")[1]],base)] = [line.strip().split(",")[10],float(line.strip().split(",")[11])]
    return return_dict

def upstream_build(fyle,ends,duplex):
    return_dict,fasta_dict = {},{}
    for record in SeqIO.parse(fyle, "fasta"):
        fasta_dict[int(record.id.split("_")[0].split(":")[1])] = str(record.seq)
    for item in ends:
        if fasta_dict[item].find(duplex[item]) != -1:
            return_dict[item] = fasta_dict[item][fasta_dict[item].find(duplex[item])-11:fasta_dict[item].find(duplex[item])]
        else:
            print 'no'
    return return_dict

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('ends',type=str,help='')
    parser.add_argument('strands',type=str,help='')
    parser.add_argument('duplex',type=str,help='')
    parser.add_argument('deltaT',type=str,help='')
    parser.add_argument('supplemental',type=str,help='')
    parser.add_argument('upstream',type=str,help='')
    parser.add_argument('outfile',type=str,help='')
    args = parser.parse_args()

    strand_list,ends,strand_dict = build_strand_list(args.strands),sorted(build_end_list(args.ends)),build_strand_dict(args.strands)
    duplex,deltaT = duplex_translator(build_duplex_dict(args.duplex),strand_list),deltaT_dict_build(args.deltaT)
    supplemental,upstream = supplemental_build(args.supplemental,strand_list),upstream_build(args.upstream,ends,duplex)
    with open(args.outfile,'w') as outp:
        outp.write('pos,strand,upstream,hairpin,hairpin_deltaG,U-tract,deltaT\n')
        for end in ends:
            outp.write(str(end)+","+str(strand_dict[end])+","+str(upstream[end])+","+str(duplex[end])+","+str(supplemental[end][1])+","+str(supplemental[end][0])+","+str(deltaT[end])+"\n")




if __name__ == '__main__':
    main()
