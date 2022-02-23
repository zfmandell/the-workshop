#!/usr/bin/env python3

import argparse
import collections

def read_ends(fyle):
	return_dict = {}
	with open(fyle,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			coord = line.strip().split(',')[0]
			strand = line.strip().split(',')[1]
			delta = float(line.strip().split(',')[5])
			info = ','.join([line.strip().split(',')[2]]+line.strip().split(',')[6:])+'\n'

			return_dict[coord+'_'+strand] = [delta,info]
	
	return return_dict


def term_crunch(one,two):
	return_list = []
	for key,value in one.items():
		delta_one = value[0]
		try:
			delta_two = two[key][0]
			return_list.append([key.split('_')[0],key.replace('_',',')+','+','.join([str(delta_one),str(delta_two),value[1]])])
		except KeyError:
			pass

	return return_list

def main():
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('condition_1',type=str,help='')
	parser.add_argument('condition_2',type=str,help='')
	parser.add_argument('output',type=str,help='')
	
	args = parser.parse_args()
	
	one = read_ends(args.condition_1)
	two = read_ends(args.condition_2)

	crunched = term_crunch(one,two)
	crunched_sorted = sorted(crunched, key=lambda x: x[0])

	A = args.condition_1.split('_')[0]
	B = args.condition_2.split('_')[0]

	with open(args.output,'w') as outp:
		outp.write('coord,strand,'+A+' d%T,'+B+' d%T,dd%T,Position,category,upstream,hairpin,downstream,dG\n')
		for item in crunched_sorted:
			outp.write(item[1])


if __name__ == '__main__':
	main()