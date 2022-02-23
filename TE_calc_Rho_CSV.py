#!/usr/bin/env python3

import argparse
import numpy as np
from scipy.signal import savgol_filter

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
				   
	return genes

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

def TE_calc(coord,cov,strand,tipe):

	if tipe == 'canonical':
	
		if strand == '+':
			up = np.median(cov[coord-100:coord])
			down = np.median(cov[coord+1:coord+500])

			if down == 0.0:
				down = 1

			if up > 10:
				TE = round(((up-down)/up)*100,2)
				if TE < 0:
					TE = 0
				return TE
			
			else:
				return 'NA'
		
		else:
			up = np.median(cov[coord-1:coord+99])
			down = np.median(cov[coord-501:coord-1])

			if down == 0.0:
				down = 1

			if up > 10:
				TE = round(((up-down)/up)*100,2)
				if TE < 0:
					TE = 0
				return TE
			
			else:
				return 'NA'

	else:
		if strand == '+':
			up = np.median(cov[coord-10:coord])
			down = np.median(cov[coord+1:coord+11])

			if down == 0.0:
				down = 1

			if up > 10:
				TE = round(((up-down)/up)*100,2)
				if TE < 0:
					TE = 0
				return TE
			
			else:
				return 'NA'
		
		else:
			up = np.median(cov[coord-1:coord+9])
			down = np.median(cov[coord-11:coord-1])

			if down == 0.0:
				down = 1

			if up > 10:
				TE = round(((up-down)/up)*100,2)
				if TE < 0:
					TE = 0
				return TE
			
			else:
				return 'NA'


def read_coverage(cov_fyle):
	#create list of coverage values, position dependent
	return_list = []
	with open(cov_fyle,"r") as inp:
		for line in inp:
			return_list.append(float(line.strip().split('\t')[2]))
	return return_list

def read_ends(fyle):
	return_dict = {}
	with open(fyle,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			coord = int(line.strip().split(',')[0])
			strand = line.strip().split(',')[1]
			return_dict[coord] = strand
	
	return return_dict

def read_terms(fyle):
	return_dict = {}
	with open(fyle,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			coord = int(line.strip().split(',')[0])
			strand = line.strip().split(',')[1]
			return_dict[coord] = strand
	
	return return_dict

def GFF3_annotate(coord,GFF3,strand):
	gene = None
	if strand == '+':
		for item in GFF3[:len(GFF3)-1]:
			if coord > item[3] and coord < GFF3[GFF3.index(item)+1][2]:
				if coord-item[3] < 1000:
					gene = item[0]+'+'+str(coord-item[3])+'nt'
				else:
					gene = 'orphan transcript'

	else:
		for item in GFF3[:len(GFF3)-1]:
			if coord < item[2] and coord > GFF3[GFF3.index(item)+1][3]:
				if item[2]-coord < 1000:
					gene = item[0]+'+'+str(item[2]-coord)+'nt'
				else:
					gene = 'orphan transcript'

	if gene == None:
		return 'orphan transcript'
	else:
		return gene


def antisense_annotate(coord,strand,GFF3):
	antisense = None
	
	if strand == '+':
		region = list(range(coord,coord+250))
	else:
		region = list(range(coord-250,coord))
	
	for item in GFF3:
		subregion = list(range(item[2],item[3]))

		if bool(set(region) & set(subregion)) == True:
			if strand != item[1]:
				antisense = 'antisense to '+item[0]
				break
			else:
				antisense = 'sense to '+item[0]
				break

	if antisense == None:
		return 'intergenic'
	else:
		return antisense

def main():
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('ends',type=str,help='')
	parser.add_argument('terms',type=str,help='')
	parser.add_argument('WT_fwd',type=str,help='')
	parser.add_argument('WT_rev',type=str,help='')
	parser.add_argument('Mutant_fwd',type=str,help='')
	parser.add_argument('Mutant_rev',type=str,help='')
	parser.add_argument('sequence',type=str,help='')
	parser.add_argument('length',type=int,help='')
	parser.add_argument('GFF3',type=str,help='')
	parser.add_argument('output',type=str,help='')
	
	args = parser.parse_args()
	
	ends = read_ends(args.ends)
	terms = read_terms(args.terms)
	WT_fwd = read_coverage(args.WT_fwd)
	WT_rev = read_coverage(args.WT_rev)
	Mutant_fwd = read_coverage(args.Mutant_fwd)
	Mutant_rev = read_coverage(args.Mutant_rev)
	genome = genome_yield(args.sequence)
	length = args.length
	GFF3 = build_GFF3(args.GFF3)
	genes_plus = [x for x in GFF3 if x[1] == '+']
	genes_minus = [x for x in GFF3 if x[1] == '-'][::-1]

	return_dict = {}
	for key,value in ends.items():
		coord = key
		strand = value
		if strand == '+':
			if coord in terms.keys():
				if terms[coord] == strand:
					category = 'intrinsic'
				else:
					category = 'canonical'
			else:
				category = 'canonical'

			if category == 'intrinsic':
				WT_TE = TE_calc(coord,WT_fwd,'+','intrinsic')
				Mutant_TE = TE_calc(coord,Mutant_fwd,'+','intrinsic')
			else:
				WT_TE = TE_calc(coord,WT_fwd,'+','canonical')
				Mutant_TE = TE_calc(coord,Mutant_fwd,'+','canonical')

			if WT_TE != 'NA' and Mutant_TE != 'NA':
				if WT_TE >= 50 and round(WT_TE-Mutant_TE,0) >= 15:
					gene = GFF3_annotate(coord,genes_plus,'+')
					antisense = antisense_annotate(coord,'+',GFF3)
					upstream = genome[coord-length:coord].replace('T','U')
					downstream = genome[coord:coord+length].replace('T','U')
					
					return_dict[coord] = [strand,gene,antisense,WT_TE,Mutant_TE,round(WT_TE-Mutant_TE,0),category,upstream,downstream]

		else:
			if coord in terms.keys():
				if terms[coord] == strand:
					category = 'intrinsic'
				else:
					category = 'canonical'
			else:
				category = 'canonical'
			
			if category == 'intrinsic':
				WT_TE = TE_calc(coord,WT_rev,'-','intrinsic')
				Mutant_TE = TE_calc(coord,Mutant_rev,'-','intrinsic')
			else:
				WT_TE = TE_calc(coord,WT_rev,'-','canonical')
				Mutant_TE = TE_calc(coord,Mutant_rev,'-','canonical')

			if WT_TE != 'NA' and Mutant_TE != 'NA':
				if WT_TE >= 50 and round(WT_TE-Mutant_TE,0) >= 15:
					gene = GFF3_annotate(coord,genes_minus,'-')
					antisense = antisense_annotate(coord,'-',GFF3)
					upstream = reverse_complement(genome[coord:coord+length]).replace('T','U')
					downstream = reverse_complement(genome[coord-length:coord]).replace('T','U')
					
					return_dict[coord] = [strand,gene,antisense,WT_TE,Mutant_TE,round(WT_TE-Mutant_TE,0),category,upstream,downstream]


	with open(args.output,'w') as outp:
		outp.write('Coord,Strand,Relative Position,Antisense to CDS,%T WT,%T Mutant,d%T,category,upstream,downstream\n')
		for key,value in return_dict.items():
			outp.write(str(key)+','+','.join(list(map(str,value)))+'\n')
	
if __name__ == '__main__':
	main()






