#!/usr/bin/env python3

import argparse
import numpy as np
import collections

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
			coord = int(line.strip().split(',')[0])
			strand = line.strip().split(',')[1]
			delta = float(line.strip().split(',')[5])
			return_dict[coord] = [strand,delta]
	
	return return_dict

def term_crunch(terms,fwd,rev):
	return_dict = {}
	for key,value in terms.items():
		if value[1] > 15.0:
			if value[0] == '+':
				if rev[key+1] != 0 and fwd[key+1] != 0:
					value.append(round(np.log2(rev[key+1]/fwd[key+1]),2))
					return_dict[key] = value
				elif rev[key+1] == 0 and fwd[key+1] != 0:
					value.append(round(np.log2(1/fwd[key+1]),2))
					return_dict[key] = value
				elif rev[key+1] != 0 and fwd[key+1] == 0:
					value.append(round(np.log2(rev[key+1]/1),2))
					return_dict[key] = value
				elif rev[key+1] == 0 and fwd[key+1] == 0:
					value.append(round(np.log2(1/1),2))
					return_dict[key] = value
			else:
				if rev[key+1] != 0 and fwd[key+1] != 0:
					value.append(round(np.log2(fwd[key+1]/rev[key+1]),2))
					return_dict[key] = value
				elif rev[key+1] != 0 and fwd[key+1] == 0:
					value.append(round(np.log2(1/rev[key+1]),2))
					return_dict[key] = value
				elif rev[key+1] == 0 and fwd[key+1] != 0:
					value.append(round(np.log2(fwd[key+1]/1),2))
					return_dict[key] = value
				elif rev[key+1] == 0 and fwd[key+1] == 0:
					value.append(round(np.log2(1/1),2))
					return_dict[key] = value
	
	return return_dict

def main():
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('terms',type=str,help='')
	parser.add_argument('fwd',type=str,help='')
	parser.add_argument('rev',type=str,help='')
	parser.add_argument('output',type=str,help='')
	
	args = parser.parse_args()
	
	terms = read_terms(args.terms)
	fwd = read_coverage(args.fwd)
	rev = read_coverage(args.rev)

	crunched = term_crunch(terms,fwd,rev)
	crunched_sorted = collections.OrderedDict(sorted(crunched.items()))

	with open(args.output,'w') as outp:
		outp.write('coord,strand,d%T,antisense\n')
		for key,value in crunched_sorted.items():
			outp.write(str(key)+','+str(value[0])+','+str(value[1])+','+str(value[2])+'\n')

if __name__ == '__main__':
	main()
