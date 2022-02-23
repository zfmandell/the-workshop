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

def build_bedgraph(fyle):
	return_dict = {}
	with open(fyle,'r') as inp:
		for line in inp:
			coord = int(line.strip().split()[2])
			delta = float(line.strip().split()[3])

			if delta > 0:
				strand = '+'
			else:
				strand = '-'
			
			return_dict[coord] = strand
	
	return return_dict

genome = genome_yield(sys.argv[1])
bed = build_bedgraph(sys.argv[2])

with open(sys.argv[3],'w') as outp:
	for key,value in bed.items():
		if value == '+':
			outp.write('>'+str(key)+'_'+str(value)+'\n'+genome[key-49:key+1].replace('T','U')+'\n')
		else:
			outp.write('>'+str(key)+'_'+str(value)+'\n'+reverse_complement(genome[key:key+50]).replace('T','U')+'\n')





