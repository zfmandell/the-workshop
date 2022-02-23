#!/usr/bin/env python3

import argparse
import numpy as np

def read_coverage(cov_fyle):
	#create list of coverage values, position dependent
	return_dict = {}
	with open(cov_fyle,"r") as inp:
		for line in inp:
			try:
				return_dict[int(line.strip().split()[0])] = float(line.strip().split()[1])
			except ValueError:
				pass
	return return_dict

def read_ends(fyle):
	return_dict = {}
	with open(fyle,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			return_dict[int(line.strip().split(',')[0])] = line.strip().split(',')[1]
	return return_dict

def ribo_crunch(strand,cov,key):
	temp = []
	if strand == '+':
		for item in range(key-25,key):
			if item in cov.keys():
				temp.append(cov[item])
	else:
		for item in range(key-1,key+25-1):
			if item in cov.keys():
				temp.append(cov[item])

	if len(temp) > 10:
		to_return = round(np.max(temp),2)
	else:
		to_return = 'NA'

	return to_return

def GFF_build(fyle,mono):
	return_dict = {}
	with open(fyle,'r') as inp:
		for line in inp:
			if len(line.strip().split('\t')) > 1:
				if line.strip().split('\t')[2] == 'CDS':
					for item in line.strip().split('\t')[8].split(';'):
						if 'gene=' in item:
							gene = item.split('=')[1]
					if gene in mono:
						strand = line.strip().split('\t')[6]
						if strand == '+':
							left = int(line.strip().split('\t')[3])
							right = int(line.strip().split('\t')[4])
						else:
							left = int(line.strip().split('\t')[3])
							right = int(line.strip().split('\t')[4])
						return_dict[gene] = [strand,left,right]

	return return_dict

def ribo_control(fwd,rev,GFF3):
	return_dict = {}
	for key,value in GFF3.items():
		if value[0] == '+':
			for item in range(value[1],value[2]):
				if item in fwd.keys():
					return_dict[item] = fwd[item]
		else:
			for item in range(value[1],value[2]):
				if item in rev.keys():
					return_dict[item] = rev[item]
	return return_dict
	
def main():
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('ribo_fwd',type=str,help='')
	parser.add_argument('ribo_rev',type=str,help='')
	parser.add_argument('ends',type=str,help='')
	parser.add_argument('GFF3',type=str,help='')
	parser.add_argument('monocistronic',type=str,help='')
	parser.add_argument('output',type=str,help='')
	
	args = parser.parse_args()

	mono = []
	with open(args.monocistronic,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			mono.append(line.strip())

	annotation = GFF_build(args.GFF3,mono)

	ends = read_ends(args.ends)
	fwd = read_coverage(args.ribo_fwd)
	rev = read_coverage(args.ribo_rev)

	final = {}
	for key,value in ends.items():
		if value == '+':
			final[key] = ribo_crunch('+',fwd,key)

		else:
			final[key] = ribo_crunch('-',rev,key)

	control = ribo_control(fwd,rev,annotation)

	with open(args.output,'w') as outp:
		outp.write('key,value\n')
		for key,value in final.items():
			outp.write(str(key)+','+str(value)+'\n')

	with open('control.csv','w') as outp:
		outp.write('key,value\n')
		for key,value in control.items():
			outp.write(str(key)+','+str(value)+'\n')

if __name__ == '__main__':
	main()
