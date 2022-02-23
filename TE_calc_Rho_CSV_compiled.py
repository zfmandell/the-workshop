#!/usr/bin/env python3

import argparse
import numpy as np
from scipy.signal import savgol_filter

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
			rel = line.strip().split(',')[2]
			anti = line.strip().split(',')[3]
			rest = ','.join(line.strip().split(',')[7:])      
			return_dict[coord] = [strand,rel,anti,rest]

	return return_dict

def main():
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('ends',type=str,help='')
	parser.add_argument('WT_fwd',type=str,help='')
	parser.add_argument('WT_rev',type=str,help='')
	parser.add_argument('Mutant_fwd',type=str,help='')
	parser.add_argument('Mutant_rev',type=str,help='')
	parser.add_argument('output',type=str,help='')
	
	args = parser.parse_args()
	
	ends = read_ends(args.ends)
	WT_fwd = read_coverage(args.WT_fwd)
	WT_rev = read_coverage(args.WT_rev)
	Mutant_fwd = read_coverage(args.Mutant_fwd)
	Mutant_rev = read_coverage(args.Mutant_rev)

	return_dict = {}
	for key,value in ends.items():
		coord = key
		strand = value[0]
		rel = value[1]
		anti = value[2]
		tipe = value[3].split(',')[0]
		rest = value[3]
		if strand == '+':
			if tipe == 'intrinsic':
				WT_TE = TE_calc(coord,WT_fwd,'+','intrinsic')
				Mutant_TE = TE_calc(coord,Mutant_fwd,'+','intrinsic')
			else:
				WT_TE = TE_calc(coord,WT_fwd,'+','canonical')
				Mutant_TE = TE_calc(coord,Mutant_fwd,'+','canonical')

			if WT_TE != 'NA' and Mutant_TE != 'NA':
				return_dict[coord] = [strand,rel,anti,WT_TE,Mutant_TE,round(WT_TE-Mutant_TE,0),rest]

		else:
			if tipe == 'intrinsic':
				WT_TE = TE_calc(coord,WT_rev,'-','intrinsic')
				Mutant_TE = TE_calc(coord,Mutant_rev,'-','intrinsic')
			else:
				WT_TE = TE_calc(coord,WT_rev,'-','canonical')
				Mutant_TE = TE_calc(coord,Mutant_rev,'-','canonical')

			if WT_TE != 'NA' and Mutant_TE != 'NA':
				return_dict[coord] = [strand,rel,anti,WT_TE,Mutant_TE,round(WT_TE-Mutant_TE,0),rest]


	with open(args.output,'w') as outp:
		outp.write('Coord,Strand,Relative Position,Antisense to CDS,%T WT,%T Mutant,d%T,category,upstream bubble density,upstream Num YC dimers,downstream bubble density,downstream Num YC dimers\n')
		for key,value in return_dict.items():
			outp.write(str(key)+','+','.join(list(map(str,value)))+'\n')
	
if __name__ == '__main__':
	main()

