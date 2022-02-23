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
		up = np.median(cov[coord-10:coord])
		down = np.median(cov[coord+1:coord+11])
		return ((up-down)/up)*100
	
	else:
		up = np.median(cov[coord+1:coord+11])
		down = np.median(cov[coord-10:coord])
		return ((up-down)/up)*100

def U_check(coord,genome,strand):
	if strand == '+':
		U = genome[coord-8:coord].replace('T','U')
	else:
		U = reverse_complement(genome[coord:coord+9]).replace('T','U')

	if (U.count('U')/len(U))*100 > 40:
		return True

def build_GFF3(fyle):
	genes,rRNA = [],[]
	with open(fyle,'r') as inp:
		for line in inp:
			if len(line.strip().split('\t')) > 1:
				if line.strip().split('\t')[2] == 'gene':
					gene = line.strip().split('\t')[8].split(';')[2].split('=')[1]
					start = int(line.strip().split('\t')[3])
					end = int(line.strip().split('\t')[4])
					strand = line.strip().split('\t')[6]
					genes.append([gene,strand,start,end])

				elif line.strip().split('\t')[2] == 'rRNA':
					start = int(line.strip().split('\t')[3])
					end = int(line.strip().split('\t')[4])
					strand = line.strip().split('\t')[6]
					rRNA.append([strand,start,end])
				   
	return genes,rRNA

"""

def build_GFF3(gff3):
	genes,rRNA = [],[]
	with open(gff3,'r') as inp:
		for line in inp:
			if len(line.strip().split('\t')) > 1:
				if line.strip().split('\t')[2] == 'CDS':
				
					CDS = line.strip().split('\t')[8].split(';')[7].split('=')[1]
					start = int(line.strip().split('\t')[3])
					end = int(line.strip().split('\t')[4])
					strand = line.strip().split('\t')[6]
					if CDS != 'ab initio prediction:AMIGene:2.0':
						genes.append([CDS,strand,start,end])
					else:
						CDS = line.strip().split('\t')[8].split(';')[6].split('=')[1]
						genes.append([CDS,strand,start,end])

				if line.strip().split('\t')[2] == 'rRNA':
					start = int(line.strip().split('\t')[3])
					end = int(line.strip().split('\t')[4])
					strand = line.strip().split('\t')[6]
					rRNA.append([strand,start,end])

	return genes,rRNA
"""

def build_csv(fyle):
	return_dict = {}
	with open(fyle,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			coord = int(line.strip().split(',')[0])
			strand = line.strip().split(',')[1]
			gene = line.strip().split(',')[2]

			return_dict[coord] = [strand,gene]

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

def term_find(csv,bedgraph):
	terms,found = {},[]
	for key,value in bedgraph.items():
		for subkey,subvalue in csv.items():
			if subkey-3 <= key <= subkey+3 and subvalue[0] == value[0]:
				terms[key] = subvalue
				found.append(key)

	return terms,found

def CDS_check(bedgraph,gff3):
	found = []
	for key,value in bedgraph.items():
		for item in gff3:
			if item[2] <= key <= item[3] and value[0] == item[1]:
				found.append(key)
	return found

def rRNA_check(bedgraph,rRNA):
	found = []
	for key,value in bedgraph.items():
		for item in rRNA:
			if item[1] <= key <= item[2] and value[0] == item[0]:
				found.append(key)
	return found

def intergenic_check(bedgraph,genes_plus,genes_minus):
	found = {}
	for key,value in bedgraph.items():
		if value[0] == '+':
			for item in genes_plus[:len(genes_plus)-1]:

				gene_up = item[0]
				start_up = int(item[2])
				stop_up = int(item[3])
				
				gene_down = genes_plus[genes_plus.index(item)+1][0]
				start_down = int(genes_plus[genes_plus.index(item)+1][2])
				stop_down = int(genes_plus[genes_plus.index(item)+1][3])

				if (key > stop_up) and (key < start_down):
					found[key] = gene_up

		else:
			for item in genes_minus[:len(genes_minus)-1]:

				gene_down = item[0]
				start_down = int(item[2])
				stop_down = int(item[3])
				
				gene_up = genes_minus[genes_minus.index(item)+1][0]
				start_up = int(genes_minus[genes_minus.index(item)+1][2])
				stop_up = int(genes_minus[genes_minus.index(item)+1][3])

				if (key < start_up) and (key > stop_down):
					found[key] = gene_up
	return found

def main():
	parser = argparse.ArgumentParser(description='takes a bedgraph file and a gff3 file, and annotates each 3-OH end within the bedgraph file')
	parser.add_argument('csv',type=str,help='csv file')
	parser.add_argument('bedgraph',type=str,help='bedgraph file')
	parser.add_argument('gff3',type=str,help='gff3 file, can be downloaded from ncbi database')
	parser.add_argument('fwd',type=str,help='fwd cov file')
	parser.add_argument('rev',type=str,help='rev cov file')
	parser.add_argument('genome',type=str,help='.fasta genome file')
   
	args = parser.parse_args()

	csv = build_csv(args.csv)
	bedgraph = build_bedgraph(args.bedgraph)
	gff3,rRNA = build_GFF3(args.gff3)
	fwd = read_coverage(args.fwd)
	rev = read_coverage(args.rev)
	genome = genome_yield(args.genome)

	genes_plus = [x for x in gff3 if x[1] == '+']
	genes_minus = [x for x in gff3 if x[1] == '-']

	terms,found = term_find(csv,bedgraph)

	for item in found:
		bedgraph.pop(item, None)

	found = CDS_check(bedgraph,gff3)

	for item in found:
		bedgraph.pop(item, None)

	found = rRNA_check(bedgraph,rRNA)

	for item in found:
		bedgraph.pop(item, None)

	"""
	genes = [x[1] for x in terms.values()]


	found = intergenic_check(bedgraph,genes_plus,genes_minus)

	for key,value in found.items():
		if value in genes:
			bedgraph.pop(key, None)

	"""

	peaks_all = [np.log10(abs(x[1])) for x in bedgraph.values()]
	perc = 10**np.percentile(peaks_all,25)

	with open('to_test.bedgraph','w') as outp:
		with open(args.bedgraph,'r') as inp:
			for line in inp:
				coord = int(line.strip().split()[2])
				delta = float(line.strip().split()[3])
				if delta > 0:
					strand = '+'
				else:
					strand = '-'

				if coord in bedgraph.keys() and abs(delta) >= perc:
					if strand == '+':
						if TE_calc(coord,fwd,'+') >= 50 and U_check(coord,genome,'+') == True:
							outp.write(line)
					else:
						if TE_calc(coord,rev,'-') >= 50 and U_check(coord,genome,'-') == True:
							outp.write(line)

	with open('matched_terms.csv','w') as outp:
		for key,value in terms.items():
			outp.write(str(key)+','+str(value)+'\n')


if __name__ == '__main__':
	main()







