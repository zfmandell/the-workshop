#!/usr/bin/env python3

import sys

def read_base(csv_file):
	return_dict = {}
	with open(csv_file,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			return_dict[int(line.strip().split(',')[0])] = line.strip().split(',')[1]
	return return_dict

def read_bedgraph(bed_fyle):
	#create list of coverage values, position dependent
	return_list = []
	with open(bed_fyle,"r") as inp:
		for line in inp:
			return_list.append(float(line.strip().split('\t')[3]))
	return return_list

base = read_base(sys.argv[1])
fwd = read_bedgraph(sys.argv[2])
rev = read_bedgraph(sys.argv[3])

final = {}
for key,value in base.items():
	if value == '+':
		final[key] = fwd[key]
	else:
		final[key] = rev[key]

with open(sys.argv[4],'w') as outp:
	outp.write('coord,Cv\n')
	for key,value in final.items():
		outp.write(str(key)+','+str(value)+'\n')




