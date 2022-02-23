#!/usr/bin/env python3

import argparse
import numpy as np

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

def TE_calc(coord,cov,strand):

	if coord < 251 or coord > len(cov) - 251:
		return 'NA'

	else:
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
	
def read_coverage(cov_fyle):
	#create list of coverage values, position dependent
	return_list = []
	with open(cov_fyle,"r") as inp:
		for line in inp:
			return_list.append(float(line.strip().split('\t')[2]))
	return return_list

def read_ends(fyle):
	return_list = []
	with open(fyle,'r') as inp:
		for line in inp:
			coord = line.strip().split()[2]
			delta = float(line.strip().split()[3])
			if delta > 0:
				strand = '+'
			else:
				strand = '-'

			return_list.append(coord+'_'+strand)
	
	return return_list

def read_terms(fyle):
	return_dict = {}
	with open(fyle,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			coord = int(line.strip().split(',')[0])
			strand = line.strip().split(',')[1]
			return_dict[coord] = strand
	
	return return_dict

def CDS_check(coord,strand,GFF3):
	for item in GFF3:
		if coord >= item[2] and coord <= item[3] and strand == item[1]:
			return True

def main():
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('ends',type=str,help='')
	parser.add_argument('terms',type=str,help='')
	parser.add_argument('WT_fwd',type=str,help='')
	parser.add_argument('WT_rev',type=str,help='')
	parser.add_argument('Mutant_fwd',type=str,help='')
	parser.add_argument('Mutant_rev',type=str,help='')
	parser.add_argument('GFF3',type=str,help='')
	parser.add_argument('output',type=str,help='')
	
	args = parser.parse_args()
	
	ends = read_ends(args.ends)
	terms = read_terms(args.terms)
	WT_fwd = read_coverage(args.WT_fwd)
	WT_rev = read_coverage(args.WT_rev)
	Mutant_fwd = read_coverage(args.Mutant_fwd)
	Mutant_rev = read_coverage(args.Mutant_rev)
	GFF3 = build_GFF3(args.GFF3)

	return_dict = {}
	for item in ends:
		coord = int(item.split('_')[0])
		strand = item.split('_')[1]
		if strand == '+':
			WT_TE = TE_calc(coord,WT_fwd,'+')
			Mutant_TE = TE_calc(coord,Mutant_fwd,'+')

			if WT_TE != 'NA' and Mutant_TE != 'NA' and CDS_check(coord,strand,GFF3) != True:
				if coord in terms.keys():
					if terms[coord] == strand:
						category = ' intrinsic term'
				else:
					category = ' non-intrinsic term'

				if WT_TE >= 25 and round(WT_TE-Mutant_TE,0) >= 25:
					return_dict[item] = [WT_TE,Mutant_TE,round(WT_TE-Mutant_TE,0),category]

		else:
			WT_TE = TE_calc(coord,WT_rev,'-')
			Mutant_TE = TE_calc(coord,Mutant_rev,'-')

			if WT_TE != 'NA' and Mutant_TE != 'NA' and CDS_check(coord,strand,GFF3) != True:
				if coord in terms.keys():
					if terms[coord] == strand:
						category = ' intrinsic term'
				else:
					category = ' non-intrinsic term'

				if WT_TE >= 25 and round(WT_TE-Mutant_TE,0) >= 25:
					return_dict[item] = [WT_TE,Mutant_TE,round(WT_TE-Mutant_TE,0),category]


	with open(args.output,'w') as outp:
		outp.write('coord,strand,WT,Mutant,Diff,category\n')
		for key,value in return_dict.items():
			coord = key.split('_')[0]
			strand = key.split('_')[1]

			outp.write(coord+','+strand+','+','.join(list(map(str,value)))+'\n')
	
if __name__ == '__main__':
	main()






