#!/usr/bin/env python3

import argparse
import statistics
import numpy as np 
import random

def term_reader(csv):
	return_dict = {}
	with open(csv,'r') as inp:
		firstline = inp.readline()
		for line in inp: 
			coord = int(line.strip().split(',')[0])
			strand = line.strip().split(',')[1]
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

def end_place(terms,nascent,strand,prune=False):
	terms,return_dict = [x for x,y in terms.items() if y == strand],{}
	if strand == '+':
		for coord in terms:
			nascent_terminal = nascent[coord-20:coord+1][::-1]
			nascent_baseline = nascent[coord-100:coord+1]

			if statistics.median(nascent_baseline) != 0 and all(v < 20 for v in nascent_terminal[-2:]) != True:
				temp = [x/statistics.median(nascent_baseline) for x in nascent_terminal]
				for index,item in enumerate(temp):
					if item == 0:
						temp[index] = 0+random.uniform(0.01,0.05)/statistics.median(nascent_baseline)
				return_dict[str(coord)+'_-'] = temp
			elif statistics.median(nascent_baseline) != 0 and all(v < 20 for v in nascent_terminal[-2:]) == True:
				if prune == True:
					pass
				else:
					temp = [x/statistics.median(nascent_baseline) for x in nascent_terminal]
					for index,item in enumerate(temp):
						if item == 0:
							temp[index] = 0+random.uniform(0.01,0.05)/statistics.median(nascent_baseline)
					return_dict[str(coord)+'_-'] = temp

	
	else:
		for coord in terms:
			nascent_terminal = nascent[coord:coord+21]
			nascent_baseline = nascent[coord:coord+101]
			
			if statistics.median(nascent_baseline) != 0 and all(v < 20 for v in nascent_terminal[-2:]) != True:
				temp = [x/statistics.median(nascent_baseline) for x in nascent_terminal]
				for index,item in enumerate(temp):
					if item == 0:
						temp[index] = 0+random.uniform(0.01,0.05)/statistics.median(nascent_baseline)
				return_dict[str(coord)+'_-'] = temp
			elif statistics.median(nascent_baseline) != 0 and all(v < 20 for v in nascent_terminal[-2:]) == True:
				if prune == True:
					pass
				else:
					temp = [x/statistics.median(nascent_baseline) for x in nascent_terminal]
					for index,item in enumerate(temp):
						if item == 0:
							temp[index] = 0+random.uniform(0.01,0.05)/statistics.median(nascent_baseline)
					return_dict[str(coord)+'_-'] = temp
		
	return return_dict

def pause_calc(RNET):
	return_list = []
	for i in range(20):
		return_list.append([x[i] for x in RNET])
	
	return [statistics.median(x) for x in return_list]

def pause_compare(WT,Mutant):
	return_list = []
	for index,value in enumerate(WT):
		return_list.append(round(np.log2(value/Mutant[index]),2))

	return return_list

def main():
	parser = argparse.ArgumentParser(description='determines 3-OH end idendity by combining Term-seq and RNET-seq information')
	parser.add_argument('Termseq',type=str,help='list of Term-seq called 3-OH ends with RNET-seq adjustments: <.csv>, no default value, must be inputted')
	parser.add_argument('RNETseq_WT',type=str,help='list of RNET-seq called 3-OH ends: <.wig>, no default value, must be inputted')
	parser.add_argument('RNETseq_Mutant',type=str,help='list of RNET-seq called 3-OH ends: <.wig>, no default value, must be inputted')
	parser.add_argument('output',type=str,help='name of output file')
   
	args = parser.parse_args()

	terms = term_reader(args.Termseq)
	
	RNET_plus_WT,RNET_minus_WT = rnet_reader(args.RNETseq_WT)
	terms_RNET_plus_WT,terms_RNET_minus_WT = end_place(terms,RNET_plus_WT,'+',True),end_place(terms,RNET_minus_WT,'-',True)
	RNET_WT = {**terms_RNET_plus_WT,**terms_RNET_minus_WT}
	RNET_plus_WT,RNET_minus_WT = [],[]
	RNET_WT_pruned = [value for key,value in RNET_WT.items()]
	final_WT = pause_calc(RNET_WT_pruned)

	RNET_plus_Mutant,RNET_minus_Mutant = rnet_reader(args.RNETseq_Mutant)
	terms_RNET_plus_Mutant,terms_RNET_minus_Mutant = end_place(terms,RNET_plus_Mutant,'+'),end_place(terms,RNET_minus_Mutant,'-')
	RNET_Mutant = {**terms_RNET_plus_Mutant,**terms_RNET_minus_Mutant}
	RNET_plus_Mutant,RNET_minus_Mutant = [],[]
	RNET_Mutant_pruned = [value for key,value in RNET_Mutant.items() if key in list(RNET_WT.keys())]
	final_Mutant = pause_calc(RNET_Mutant_pruned)

	print(len(RNET_WT_pruned),len(RNET_Mutant_pruned))
	final = pause_compare(final_WT,final_Mutant)

	with open(args.output,'w') as outp:
		outp.write('Pos,Value\n')
		q=1
		for item in final:
			outp.write(str(q)+','+str(item)+'\n')
			q+=1
	
if __name__ == '__main__':
	main()