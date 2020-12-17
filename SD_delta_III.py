#!/usr/bin/env python

from __future__ import division
from itertools import islice
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import argparse


def read_fasta(genome_fasta):
    fasta_sequences,fasta_dict =SeqIO.parse(open(genome_fasta),'fasta'),{}
    for fasta in fasta_sequences:
        fasta_dict[str(fasta.id)] = str(fasta.seq)
    return fasta_dict

def filter_dictonary(in_dict,filter_dict):
    '''Removes entries from a dictionary based on coverage'''
    for k in in_dict.keys():
        if k not in filter_dict:
            del in_dict[k]

def get_covered_transcripts(coverage_fyle):
    '''Reads in a non standard overlapped coverage file from SD_delta'''
    temp = []
    info = {}
    with open(coverage_fyle,'r') as f:
        for line in f:
            line = line.strip()
            line = line.split(",")
            temp.append(line)

        strongest = max(temp,key=len)

        return_temp  = [x.strip(' ') for x in strongest]
        return_temp  = [x.strip("'") for x in return_temp]

        for item in return_temp:
            info[item] = None
    return info

def mod_dict(reactvity_dict):
    'outputs -20<x<+10 window for all possible transcripts'

    output_dict = {}

    for key, value in reactvity_dict.iteritems():
        contents = key.split("|")

        new_contents = []

        for item in contents:
            if ":" in str(item):
                new_contents.append(item.split(":"))
            else:
                new_contents.append(item)


        if str(new_contents[1][0]) == 'FPUTR' and int(new_contents[1][1]) == 20:

            output_dict[key] = value[0:31]

        elif str(new_contents[1][0]) == 'FPUTR' and int(new_contents[1][1]) > 20:
            distance = int(new_contents[1][1])-20

            output_dict[key] = value[distance:distance+31]

    return output_dict

def averaged_fasta(fasta_dict,name):

    with open(name,'w') as f:
        f.write("Pos,A,T,C,G\n")

        q = 0
        i =-20

        values = fasta_dict.values()
        length = len(values[0])

        while q < length:
            temp  = [item[q] for item in values]

            width = len(temp)

            A = temp.count("A")
            T = temp.count("T")
            G = temp.count("G")
            C = temp.count("C")

            A_perc = int(A)/float(width)
            T_perc = int(T)/float(width)
            G_perc = int(G)/float(width)
            C_perc = int(C)/float(width)

            f.write(str(i)+","+str(A_perc)+","+str(T_perc)+","+str(C_perc)+","+str(G_perc))
            f.write("\n")

            q+= 1
            i+= 1


def main():
    parser = argparse.ArgumentParser(description='processing from SD_delta.py')
    parser.add_argument('restrict',default = None, help = '<.txt > Limit analysis to these specific transcripts')
    parser.add_argument('fasta',type=str, default=None, help = '<.txt > desired fasta')
    parser.add_argument('-name',type=str, default='test.txt', help = '<.txt > desired name')

    args = parser.parse_args()

    fasta = read_fasta(args.fasta)

    covered = get_covered_transcripts(args.restrict)
    filter_dictonary(fasta,covered)

    SD_fasta = mod_dict(fasta)

    averaged_name = str(args.name)+".csv"

    averaged_fasta(SD_fasta,averaged_name)


if __name__ == '__main__':
    main()
