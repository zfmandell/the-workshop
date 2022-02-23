#!/usr/bin/env python3

import sys

def genome_yield(fasta_name):
    seq = ''
    with open(fasta_name) as inp:
        for line in inp:
            if line[0] != '>':
                seq = seq+str(line.strip())
    return seq

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases

def coord_return(fyle):
    return_dict = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            return_dict[int(line.strip().split(',')[0])] = line.strip().split(',')[1]
    return return_dict

genome = genome_yield(sys.argv[1])
coords = coord_return(sys.argv[2])
length = int(sys.argv[3])

with open('sequences.fa','w') as outp:
    for key,value in coords.items():
        if value == '+':
            outp.write('>'+str(key)+'_'+str(value)+'\n'+genome[key-length:key].replace('T','U')+'\n')
        else:
             outp.write('>'+str(key)+'_'+str(value)+'\n'+reverse_complement(genome[key-1:key+length-1]).replace('T','U')+'\n')

