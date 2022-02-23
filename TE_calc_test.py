#!/usr/bin/env python3

import argparse
import numpy as np
from operator import itemgetter

def TE_calc(coord,cov,strand,tipe):
	if tipe == 'term':
		if strand == '+':
			up = np.median(cov[coord-10:coord])
			down = np.median(cov[coord+1:coord+10])

			if up < 10:
				TE = 'NA'
			else:	
				TE = round(((up-down)/up)*100,2)
				if TE < 0.0:
					TE = 0.0 
			return [str(up),str(TE)]
		
		else:
			up = np.median(cov[coord-1:coord+9])
			down = np.median(cov[coord-11:coord-1])

			if up < 10:
				TE = 'NA'
			else:
				TE = round(((up-down)/up)*100,2)
				if TE < 0.0:
					TE = 0.0 
			return [str(up),str(TE)]
	
	elif tipe == 'rnet':
		if strand == '+':
			up = np.median(cov[coord-20:coord])
			down = np.median(cov[coord+9:coord+29])

			if up < 10:
				TE = 'NA'
			else:

				if down == 0:
					down = 1
				
				TE = round((up/down),2)
				
				if TE < 0.0:
					TE = 0.0 
			return [str(up),str(TE)]
		
		else:
			up = np.median(cov[coord-1:coord+19])
			down = np.median(cov[coord-30:coord-10])

			if up < 10:
				TE = 'NA'
			else:

				if down == 0:
					down = 1

				TE = round((up/down),2)
				
				if TE < 0.0:
					TE = 0.0 
			return [str(up),str(TE)]
	else:
		print('choose term or rnet')

def read_coverage(cov_fyle):
	#create list of coverage values, position dependent
	return_list = []
	with open(cov_fyle,"r") as inp:
		for line in inp:
			return_list.append(float(line.strip().split('\t')[2]))
	return return_list

def read_terms(fyle):
	return_list = []
	with open(fyle,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			coord = line.strip().split(',')[0]
			strand = line.strip().split(',')[1]

			return_list.append(coord+'_'+strand)
	
	return return_list

def main():
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('terms',type=str,help='')
	parser.add_argument('fwd',type=str,help='')
	parser.add_argument('rev',type=str,help='')
	parser.add_argument('type',type=str,help='')
	parser.add_argument('output',type=str,help='')

	args = parser.parse_args()

	ends = read_terms(args.terms)
	fwd = read_coverage(args.fwd)
	rev = read_coverage(args.rev)

	final_list = []
	for item in ends:
		key = int(item.split('_')[0])
		value = item.split('_')[1]
		if value == '+':
			final_list.append([key,value,TE_calc(key,fwd,'+',args.type)])
		else:
			final_list.append([key,value,TE_calc(key,rev,'-',args.type)])

	final_list_sorted = sorted(final_list, key=itemgetter(0))

	with open(args.output,'w') as outp:
		outp.write('POT,Strand,up,TE\n')
		for item in final_list_sorted:
			coord = str(item[0])
			strand = item[1]
			up = item[2][0]
			TE = item[2][1]
			outp.write(coord+','+strand+','+up+','+TE+'\n')
		
if __name__ == '__main__':
	main()