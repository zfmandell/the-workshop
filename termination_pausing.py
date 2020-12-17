#!/usr/bin/env python3

import argparse
import statistics
import numpy as np 
import random

def restricter(csv):
	return_list = []
	with open(csv,'r') as inp:
		firstline = inp.readline()
		for line in inp: 
			coord = line.strip().split(',')[0]
			strand = line.strip().split(',')[1]
			return_list.append(coord+'_'+strand)

	return return_list

def term_reader(csv,restrictor):
	return_dict = {}
	with open(csv,'r') as inp:
		firstline = inp.readline()
		for line in inp: 
			coord = line.strip().split(',')[0]
			strand = line.strip().split(',')[1]
			if coord+'_'+strand in restrictor:
				coord = int(coord)
				dist = int(line.strip().split(',')[2])
				if strand == '+':
					return_dict[coord+dist] = strand
				else:
					return_dict[coord-dist] = strand
	
	return return_dict

def rnet_reader(bw):
	return_plus,return_minus = [],[]
	with open(bw,'r') as inp:
		for line in inp:
			count = int(line.strip().split('\t')[1])
			if count > 0:
				return_plus.append(count)
				return_minus.append(0)
			if count < 0:
				return_plus.append(0)
				return_minus.append(abs(count))
			elif count == 0:
				return_plus.append(0)
				return_minus.append(0)

	to_add = [0] * 101
	return_plus = to_add + return_plus + to_add
	return_minus = to_add + return_minus + to_add

	return return_plus,return_minus

def end_place(released,nascent,strand):
	ends,return_dict = [x for x,y in released.items() if y == strand],{}
	if strand == '+':
		for coord in ends:
			nascent_terminal = nascent[coord]
			nascent_baseline = nascent[coord-100:coord+1]

			if nascent_terminal == 0:
				nascent_terminal = 0.01
			
			if statistics.median(nascent_baseline) == 0 :
				pass
			else:
				nascent_normal = nascent_terminal/statistics.median(nascent_baseline)
				return_dict[str(coord)+'_+'] = nascent_normal
				
	else:
		for coord in ends:
			nascent_terminal = nascent[coord]
			nascent_baseline = nascent[coord:coord+101]

			if nascent_terminal == 0:
				nascent_terminal = 0.01

			if statistics.median(nascent_baseline) == 0 :
				pass
			else:
				nascent_normal = nascent_terminal/statistics.median(nascent_baseline)
				return_dict[str(coord)+'_+'] = nascent_normal
	
	return return_dict

def pause_calc(WT,Mutant):
	return_dict = {}
	for key,value in WT.items():
		temp = []
		try:
			return_dict[key] = round(np.log2(value/Mutant[key]),2)
		except KeyError:
			pass

	return return_dict

def main():
	parser = argparse.ArgumentParser(description='determines 3-OH end idendity by combining Term-seq and RNET-seq information')
	parser.add_argument('category',type=str,help='list of Term-seq called 3-OH ends of a particular category: <.csv>, no default value, must be inputted')
	parser.add_argument('ends',type=str,help='total list of Term-seq called 3-OH ends with RNET-seq adjustments: <.csv>, no default value, must be inputted')
	parser.add_argument('RNETseq_WT',type=str,help='list of RNET-seq called 3-OH ends in WT: <.wig>, no default value, must be inputted')
	parser.add_argument('RNETseq_Mutant',type=str,help='list of RNET-seq called 3-OH ends in Mutant: <.wig>, no default value, must be inputted')
	parser.add_argument('output',type=str,help='name of output file')
   
	args = parser.parse_args()

	categorical = restricter(args.category)
	print(len(categorical))
	terms = term_reader(args.ends,categorical)
	print(len(terms.keys()))
	
	WT_plus,WT_minus = rnet_reader(args.RNETseq_WT)
	terms_WT_plus,terms_WT_minus = end_place(terms,WT_plus,'+'),end_place(terms,WT_minus,'-')
	WT = {**terms_WT_plus,**terms_WT_minus}
	WT_plus,WT_minus = [],[]
	

	Mutant_plus,Mutant_minus = rnet_reader(args.RNETseq_Mutant)
	terms_Mutant_plus,terms_Mutant_minus = end_place(terms,Mutant_plus,'+'),end_place(terms,Mutant_minus,'-')
	Mutant = {**terms_Mutant_plus,**terms_Mutant_minus}
	Mutant_plus,Mutant_minus = [],[]

	final = pause_calc(WT,Mutant)

	with open(args.output,'w') as outp:
		outp.write('coord,strand,log2\n')
		for key,value in final.items():
			outp.write(str(key.split('_')[0])+','+str(key.split('_')[1])+','+str(value)+'\n')
	
if __name__ == '__main__':
	main()