#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype
from plotnine import *

def GFF_build(fyle):
	return_list = []
	with open(fyle,'r') as inp:
		for line in inp:
			if len(line.strip().split('\t')) > 1:
				if line.strip().split('\t')[2] == 'gene':
					start = int(line.strip().split('\t')[3])
					end = int(line.strip().split('\t')[4])
					strand = line.strip().split('\t')[6]
					return_list.append([start,end,strand])
	
	return return_list

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
		firstline = inp.readline()
		for line in inp:
			coord = line.strip().split(',')[0]
			strand = line.strip().split(',')[1]
			return_list.append(coord+'_'+strand)
	
	return return_list

def TPM_scaling(fwd,rev,genes):
	RPK = []
	for value in genes:
		if value[2] == '+':
			RPK.append(sum(fwd[value[0]:value[1]+1])/(abs(value[0]-value[1])+1/1000))
		else:
			RPK.append(sum(rev[value[0]:value[1]+1])/(abs(value[0]-value[1])+1/1000))

	return sum(RPK)/1000000

def TPM_crunch(ends,upstream,downstream,window,step,fwd,rev,scaling):
	final_list = []
	for item in ends:
		coord = int(item.split('_')[0])
		strand = item.split('_')[1]
		if coord-upstream > 1 and coord+downstream < len(fwd): 
			TPM = []
			if strand == '+':
				bins = [fwd[i:i+window] for i in list(range(coord-upstream,coord+downstream+1,step))]
				for subitem in bins:
					TPM.append(round((sum(subitem)/(window/1000))/scaling,1))
			else:
				bins = [rev[i:i+window] for i in list(range(coord-downstream,coord+upstream+1,step))][::-1]
				for subitem in bins:
					TPM.append(round((sum(subitem)/(window/1000))/scaling,1))
			final_list.append(TPM)

	median_list = []
	for item in range(0,len(final_list[0])):
		median_list.append(round(np.log10(np.median([x[item] for x in final_list])),1))

	return median_list

def main():
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('ends',type=str,help='')
	parser.add_argument('GFF3',type=str,help='')
	parser.add_argument('fwd_A',type=str,help='')
	parser.add_argument('rev_A',type=str,help='')
	parser.add_argument('fwd_B',type=str,help='')
	parser.add_argument('rev_B',type=str,help='')
	parser.add_argument('condition_A',type=str,help='')
	parser.add_argument('condition_B',type=str,help='')
	parser.add_argument('upstream',type=int,help='')
	parser.add_argument('downstream',type=int,help='')
	parser.add_argument('window',type=int,help='')
	parser.add_argument('step',type=int,help='')
	parser.add_argument('output',type=str,help='')

	args = parser.parse_args()

	ends = read_ends(args.ends)
	GFF3 = GFF_build(args.GFF3)
	fwd_A = read_coverage(args.fwd_A)
	rev_A = read_coverage(args.rev_A)
	fwd_B = read_coverage(args.fwd_B)
	rev_B = read_coverage(args.rev_B)
	scaling_A = TPM_scaling(fwd_A,rev_A,GFF3)
	scaling_B = TPM_scaling(fwd_B,rev_B,GFF3)

	median_list_A = TPM_crunch(ends,args.upstream,args.downstream,args.window,args.step,fwd_A,rev_A,scaling_A)
	median_list_B = TPM_crunch(ends,args.upstream,args.downstream,args.window,args.step,fwd_B,rev_B,scaling_B)

	all_loc = []
	current = args.upstream*-1
	all_loc.append(current)
	while current < args.downstream:
		current += args.step
		all_loc.append(current)

	all_cond = [args.condition_A]*len(all_loc)+[args.condition_B]*len(all_loc)
	all_loc = all_loc*2
	all_median = median_list_A+median_list_B

	final_for_df = {'distance':all_loc,'condition':all_cond,'TPM':all_median}

	df = pd.DataFrame(final_for_df)
	
	plot = (ggplot(df)
	 	+ aes(x='distance',y='TPM',group='condition',color='condition')
	 	+ scale_color_grey() 
	 	+ theme_classic()
	 	+ geom_line(linetype="dashed", size=1.2)
	 	+ geom_point(color="red", size=3)        
	 	+ labs(x='distance', y='log10 TPM')
	 	)
	plot.save(args.output, height=6, width=8)
 	
	

if __name__ == '__main__':
	main()