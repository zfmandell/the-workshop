#!/usr/bin/env python3

import argparse
import numpy as np
import operator

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

def read_coverage(cov_fyle):
	#create list of coverage values, position dependent
	return_list = []
	with open(cov_fyle,"r") as inp:
		for line in inp:
			return_list.append(float(line.strip().split('\t')[2]))
	return return_list

def TE_calc(coord,cov,strand):
	if strand == '+':
		up = np.median(cov[coord-20:coord])
		down = np.median(cov[coord+1:coord+21])
		return ((up-down)/up)*100
	
	else:
		up = np.median(cov[coord+1:coord+21])
		down = np.median(cov[coord-20:coord])
		return ((up-down)/up)*100

def U_check(coord,genome,strand):
	if strand == '+':
		U = genome[coord-8:coord].replace('T','U')
	else:
		U = reverse_complement(genome[coord:coord+9]).replace('T','U')

	if (U.count('U')/len(U))*100 > 40:
		return True


def build_csv(fyle):
	return_dict = {}
	with open(fyle,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			coord = int(line.strip().split(',')[0])
			strand = line.strip().split(',')[1]
			gene = line.strip().split(',')[2]

			return_dict[coord] = strand

	return return_dict


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
			

			return_dict[coord] = [strand,delta]
	
	return return_dict

def rRNA_check(bedgraph,rRNA):
	found = []
	for key,value in bedgraph.items():
		for item in rRNA:
			if item[1] <= key <= item[2] and value[0] == item[0]:
				found.append(key)
	return found

def build_GFF3(gff3):
	rRNA = []
	with open(gff3,'r') as inp:
		for line in inp:
			if len(line.strip().split('\t')) > 1:
				if line.strip().split('\t')[2] == 'rRNA':
					start = int(line.strip().split('\t')[3])
					end = int(line.strip().split('\t')[4])
					strand = line.strip().split('\t')[6]
					rRNA.append([strand,start,end])

	return rRNA

def cov_check(coord,cov,strand):
	if strand == '+':
		up = np.median(cov[coord-20:coord])
	else:
		up = np.median(cov[coord+1:coord+21])

	if up >= 10:
		return True

def main():
	parser = argparse.ArgumentParser(description='takes a bedgraph file and a gff3 file, and annotates each 3-OH end within the bedgraph file')
	parser.add_argument('bedgraph',type=str,help='bedgraph file')
	parser.add_argument('gff3',type=str,help='gff3 file, can be downloaded from ncbi database')
	parser.add_argument('fwd_RNA',type=str,help='fwd cov file')
	parser.add_argument('rev_RNA',type=str,help='rev cov file')
	parser.add_argument('genome',type=str,help='.fasta genome file')
   
	args = parser.parse_args()

	bedgraph = build_bedgraph(args.bedgraph)
	rRNA = build_GFF3(args.gff3)
	fwd = read_coverage(args.fwd_RNA)
	rev = read_coverage(args.rev_RNA)

	genome = genome_yield(args.genome)

	found = rRNA_check(bedgraph,rRNA)

	for item in found:
		bedgraph.pop(item, None)

	with open('to_test_RNET.bedgraph','w') as outp:
		with open(args.bedgraph,'r') as inp:
			for line in inp:
				coord = int(line.strip().split()[2])
				delta = float(line.strip().split()[3])
				if delta > 0:
					strand = '+'
				else:
					strand = '-'

				if coord in bedgraph.keys():
					if strand == '+':
						if TE_calc(coord,fwd,'+') >= 50 and U_check(coord,genome,'+') == True and cov_check(coord,fwd,'+') == True:
							outp.write(line)
					else:
						if TE_calc(coord,rev,'-') >= 50 and U_check(coord,genome,'-') == True and cov_check(coord,rev,'-') == True:
							outp.write(line)

if __name__ == '__main__':
	main()


