#!/usr/bin/env python3

import argparse
import numpy as np

def read_wig(fyle):
	return_list = []
	with open(fyle,'r') as inp:
		for line in inp:
			return_list.append(int(line.strip().split()[1]))
	return return_list

def normalize_wig(wig):
	positive = [x if x > 0 else 0 for x in wig]
	negative = [x if x < 0 else 0 for x in wig]

	p_normal,n_normal = [],[]
	for index,item in enumerate(positive):
		upstream = positive[index-250:index]
		if len(upstream) == 250:
			if np.median(upstream) != 0:
				p_normal.append(item/np.median(upstream))
			else:
				p_normal.append(item)
		else:
			p_normal.append(0)
	for index,item in enumerate(negative):
		upstream = negative[index:index+250]
		if len(upstream) == 250:
			if np.median(upstream) != 0:
				n_normal.append(item/np.median(upstream))
			else:
				n_normal.append(item)
		else:
			n_normal.append(0)

	return p_normal,n_normal

def main():
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('wig',type=str,help='')
	parser.add_argument('output',type=str,help='')
	
	args = parser.parse_args()

	wig = read_wig(args.wig)
	p,n = normalize_wig(wig)

	with open(args.output,'w') as outp:
		q = 100
		for item in p:
			outp.write('NC_000964.3\t'+str(q)+'\t'+str(q+1)+'\t'+str(item)+'\n')
			q+=1
		q = 100
		for item in n:
			outp.write('NC_000964.3\t'+str(q)+'\t'+str(q+1)+'\t'+str(-item)+'\n')
			q+=1

if __name__ == '__main__':
	main()