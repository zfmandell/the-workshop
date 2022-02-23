#!/usr/bin/env python3

import argparse
import numpy as np
from operator import itemgetter
import codecs

def TE_calc(coord,cov,strand):
	
	if strand == '+':
		up = np.median(cov[coord-10:coord])
		down = np.median(cov[coord+1:coord+10])

		if up >= 10:
			TE = round(((up-down)/up)*100,2)
			if TE < 0.0:
				TE = 0.0 
			return TE
		else:
			return None

	else:
		up = np.median(cov[coord-1:coord+9])
		down = np.median(cov[coord-11:coord-1])

		if up >= 10:
			TE = round(((up-down)/up)*100,2)
			if TE < 0.0:
				TE = 0.0 
			return TE
		else:
			return None
	
def read_coverage(cov_fyle):
	#create list of coverage values, position dependent
	return_list = []
	with open(cov_fyle,"r") as inp:
		for line in inp:
			return_list.append(float(line.strip().split('\t')[2]))
	return return_list

def read_terms(fyle):
	return_dict = {}
	with open(fyle,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			coord = line.strip().split(',')[0]
			strand = line.strip().split(',')[1]
			pos = line.strip().split(',')[2]
			A = line.strip().split(',')[3]
			hairpin = line.strip().split(',')[4]
			U = line.strip().split(',')[5]
			G  = line.strip().split(',')[6]
			return_dict[coord+'_'+strand] = [pos,A,hairpin,U,G]
	
	return return_dict

def main():
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('terms',type=str,help='')
	parser.add_argument('WT_fwd',type=str,help='')
	parser.add_argument('WT_rev',type=str,help='')
	parser.add_argument('Mutant_fwd',type=str,help='')
	parser.add_argument('Mutant_rev',type=str,help='')
	parser.add_argument('output',type=str,help='')

	args = parser.parse_args()

	ends = read_terms(args.terms)
	WT_fwd = read_coverage(args.WT_fwd)
	WT_rev = read_coverage(args.WT_rev)
	Mutant_fwd = read_coverage(args.Mutant_fwd)
	Mutant_rev = read_coverage(args.Mutant_rev)

	final_list = []
	for item in ends.keys():
		key = int(item.split('_')[0])
		value = item.split('_')[1]
		if value == '+':
			TE_WT = TE_calc(key,WT_fwd,'+')
			TE_Mutant = TE_calc(key,Mutant_fwd,'+')
			
			if TE_WT != None and TE_Mutant != None:
				if TE_WT >= 50.0:
					final_list.append([key,value,ends[item][0],TE_WT,TE_Mutant,round((TE_WT-TE_Mutant),0),ends[item][1],ends[item][2],ends[item][3],ends[item][4]])
		else:
			TE_WT = TE_calc(key,WT_rev,'-')
			TE_Mutant = TE_calc(key,Mutant_rev,'-')
			
			if TE_WT != None and TE_Mutant != None:
				if TE_WT >= 50.0:
					final_list.append([key,value,ends[item][0],TE_WT,TE_Mutant,round((TE_WT-TE_Mutant),0),ends[item][1],ends[item][2],ends[item][3],ends[item][4]])

	final_list_sorted = sorted(final_list, key=itemgetter(0))

	with open(args.output,'w',encoding='utf-8',newline='') as outp:
		outp.write('POT,Strand,Relative Position,%T WT,%T Mutant, d%T, Upstream Sequence, Predicted Hairpin, Downstream Sequence, dG-Hairpin (kcal/mol)\n')
		for item in final_list_sorted:
			outp.write(','.join(list(map(str,item)))+'\n')

		
if __name__ == '__main__':
	main()




	