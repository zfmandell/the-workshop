#!/usr/bin/env python3

import argparse

def read_new(fyle):
	return_list = []
	with open(fyle,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			coord = line.strip().split()[0].split(',')[0]
			strand = line.strip().split()[0].split(',')[1].strip('[')[1]
			return_list.append(coord+'_'+strand)
	return return_list

def read_old(fyle):
	return_dict = {}
	with open(fyle,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			coord = line.strip().split(',')[0]
			strand = line.strip().split(',')[1]
			info = ','.join(line.strip().split(',')[2:])
			return_dict[coord+'_'+strand] = info
	return return_dict,firstline 

def main():
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('new',type=str,help='')
	parser.add_argument('old',type=str,help='')
	parser.add_argument('output',type=str,help='')
   
	args = parser.parse_args()

	old,first = read_old(args.old)
	new = read_new(args.new)

	final = {}
	for item in new:
		coord = int(item.split('_')[0])
		strand = item.split('_')[1]

		for key,value in old.items():
			subcoord = int(key.split('_')[0])
			substrand = key.split('_')[1]

			if (subcoord-3 <= coord <= subcoord+3) and (strand == substrand):
				final[item] = value

	with open(args.output,'w') as outp:
		outp.write(first)
		for key,value in final.items():
			coord = key.split('_')[0]
			strand = key.split('_')[1]
			outp.write(coord+','+strand+','+value+'\n')

if __name__ == '__main__':
	main()







