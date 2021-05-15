#!/usr/bin/env python3

import argparse
import numpy as np
import csv
from collections import OrderedDict

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
							left = int(line.strip().split('\t')[3])-50
							right = int(line.strip().split('\t')[4])
						else:
							left = int(line.strip().split('\t')[3])
							right = int(line.strip().split('\t')[4])+50
						return_dict[gene] = [strand,left,right]

	return return_dict

def read_bedgraph(bed_fyle):
	#create list of coverage values, position dependent
	return_dict = {}
	with open(bed_fyle,"r") as inp:
		for line in inp:
			coord = int(line.strip().split('\t')[2])
			delta = float(line.strip().split('\t')[3])
			if delta == 0.0:
				delta = 0.01
			return_dict[coord] = delta
	return return_dict

def bed_compare(base,query,annotation,strand):
	return_dict = {}
	for key,value in annotation.items():
		if value[0] == strand:
			for item in range(value[1],value[2]+1):
				FC = round(np.log2(query[item]/base[item]),2)
				if FC > 0:
					return_dict[item] = FC
				else:
					return_dict[item] = 0
	return return_dict

def main():
	parser = argparse.ArgumentParser(description='takes two pairs of delta files, outputs set of non-zero log2FC across all transcripts. lengths equalized with NA')
	parser.add_argument('GFF3',type=str,help='GFF3 annotation file')
	parser.add_argument('monocistronic',type=str,help='list of monocistronic genes')
	parser.add_argument('base_fwd',type=str,help='fwd strand base fie')
	parser.add_argument('base_rev',type=str,help='rev strand base fie')
	parser.add_argument('query_fwd',type=str,help='fwd strand query fie')
	parser.add_argument('query_rev',type=str,help='rev strand query fie')
	parser.add_argument('output',type=str,help='base name out output files')

	args = parser.parse_args()

	mono = []
	with open(args.monocistronic,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			mono.append(line.strip())

	annotation = GFF_build(args.GFF3,mono)

	base_fwd = read_bedgraph(args.base_fwd)
	base_rev = read_bedgraph(args.base_rev)
	query_fwd = read_bedgraph(args.query_fwd)
	query_rev = read_bedgraph(args.query_rev)

	compare_fwd = bed_compare(base_fwd,query_fwd,annotation,'+')
	compare_rev = bed_compare(base_rev,query_rev,annotation,'-')
	compared = OrderedDict(sorted({**compare_fwd , **compare_rev }.items()))

	with open(args.output+'.bedgraph','w') as outp:
		for key,value in compared.items():
			outp.write('NC_000964.3\t'+str(key-1)+'\t'+str(key)+'\t'+str(value)+'\n')

	transcripts = {}
	for key,value in annotation.items():
		temp = []
		for item in range(value[1],value[2]+1):
			temp.append(compared[item])
		if value[0] == '+':
			thresh = [x for x in temp if abs(x) > 2]
			transcripts[key] = thresh
		else:
			thresh = [x for x in temp[::-1] if abs(x) > 2]
			transcripts[key] = thresh

	lengths = {}
	for key,value in transcripts.items():
		try:
			lengths[key] = sum(value)/len(value)
		except ZeroDivisionError:
			pass
	sorted_lengths = {k: v for k, v in sorted(lengths.items(), key=lambda item: item[1],reverse=True)}
	

	for key in list(sorted_lengths.keys())[:25]:
		print(key)


if __name__ == '__main__':
	main()






