#!/usr/bin/env python3

import argparse

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

def GFF_build(fyle):
	return_dict = {}
	with open(fyle,'r') as inp:
		for line in inp:
			if len(line.strip().split('\t')) > 1:
				if line.strip().split('\t')[2] == 'CDS':
					for item in line.strip().split('\t')[8].split(';'):
						if 'gene=' in item:
							gene = item.split('=')[1]
							strand = line.strip().split('\t')[6]
							left = int(line.strip().split('\t')[3])
							right = int(line.strip().split('\t')[4])
							return_dict[gene] = [strand,left,right]

	return return_dict

def seq_yield(genome,left,right,strand):
	if strand == '+':
		return genome[left:right+1]
	else:
		return reverse_complement(genome[left:right+1])

def main():
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('GFF3',type=str,help='GFF3 annotation file')
	parser.add_argument('fasta',type=str,help='list of monocistronic genes')
	parser.add_argument('output',type=str,help='name out output file')

	args = parser.parse_args()

	fasta = genome_yield(args.fasta)
	annotation = GFF_build(args.GFF3)

	final = {}
	for key,value in annotation.items():
		final['>'+str(key)] = seq_yield(fasta,value[1],value[2],value[0])

	with open(args.output,'w') as outp:
		for key,value in final.items():
			outp.write(key+'\n'+value+'\n')

if __name__ == '__main__':
	main()