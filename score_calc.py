#!/usr/bin/env python3

import argparse
import numpy as np
from operator import itemgetter

def read_terms(fyle):
	return_list = []
	with open(fyle,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			coord = line.strip().split(',')[0]
			strand = line.strip().split(',')[1]

			return_list.append(coord+'_'+strand)
	
	return return_list

def read_wig(fyle):
	return_list = []
	with open(fyle,'r') as inp:
		for line in inp:
			return_list.append(int(line.strip().split()[1]))
	return return_list

def score_calc(coord,strand,wig):
	to_add = [-5,-4,-3,-2,-1,0,1,2,3,4,5]
	coord = coord-100
	end_window = wig[coord-6:coord+5]

	if strand == '+':
		end_plus = [n if n > 0 else 0 for n in end_window]
		ind = end_plus.index(max(end_plus))
		coord = coord+to_add[ind]

		large_window = wig[coord-250:coord]
		large_plus = [abs(n) if n > 0 else 0 for n in large_window]
		
		Sum = np.log10(np.sum(large_plus))

		if Sum == 0.0:
			Sum = 1
		
		return round((large_plus[-1]/Sum),2)
	
	else:
		end_minus = [n if n < 0 else 0 for n in end_window]
		ind = end_minus.index(min(end_minus))
		coord = coord+to_add[ind]

		large_window = wig[coord-1:coord+249]
		large_minus = [abs(n) if n < 0 else 0 for n in large_window]

		Sum = np.log10(np.sum(large_minus))

		if Sum == 0.0:
			Sum = 1

		return round((large_minus[-1]/Sum),2)


def main():
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('terms',type=str,help='')
	parser.add_argument('wt',type=str,help='')
	parser.add_argument('mutant',type=str,help='')
	parser.add_argument('output',type=str,help='')
	
	args = parser.parse_args()

	terms = read_terms(args.terms)
	wt = read_wig(args.wt)
	mutant = read_wig(args.mutant)

	final_list = []
	for item in terms:
		coord = int(item.split('_')[0])
		strand = item.split('_')[1]


		wt_score = score_calc(coord,strand,wt)
		mutant_score = score_calc(coord,strand,mutant)

		if wt_score > 0:
			if mutant_score == 0.0:
				mutant_score = 1
			final_list.append([coord,strand,wt_score,mutant_score,round((wt_score/mutant_score),2)])
		else:
			final_list.append([coord,strand,'NA','NA','NA'])
	
	final_list_sorted = sorted(final_list, key=itemgetter(0))

	with open(args.output,'w') as outp:
		outp.write('POT,strand,wt,mutant,fc\n')
		for item in final_list_sorted:
			if item[0] == 4179160:
				print('yybG '+str(item[4]))
			elif item[0] == 592246:
				print('nap '+str(item[4]))
			elif item[0] == 3183216:
				print('yuae '+str(item[4]))
			coord = str(item[0])
			strand = item[1]
			wt = str(item[2])
			mutant = str(item[3])
			fc = str(item[4])
			outp.write(coord+','+strand+','+wt+','+mutant+','+fc+'\n')
	
		
if __name__ == '__main__':
	main()

