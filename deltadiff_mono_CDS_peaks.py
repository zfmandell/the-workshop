#!/usr/bin/env python3

import argparse
import numpy as np
from scipy.signal import argrelextrema
from glob import glob

def read_bedgraph(bed_fyle):
	#create list of coverage values, position dependent
	return_list = []
	with open(bed_fyle,"r") as inp:
		for line in inp:
			return_list.append(float(line.strip().split('\t')[3]))
	return return_list

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

def compare(base,query):
	diffs = np.asarray([a / b for a, b in zip(base, query)])
	minima = argrelextrema(diffs, np.less)[0].tolist()
	max_min = sorted(minima)
	
	return [enumerate(diffs),max_min]

def CDS_thresh(lst,GFF3,strand):
	return_list = []
	for item in lst:
		if queried(item+1,GFF3,strand) == True:
			return_list.append(item)

	return return_list

def queried(number,GFF3,strand):
	for key,value in GFF3.items():
		if (strand == value[0]) and (value[1] <= number <= value[2]):
			return True

def main():
	parser = argparse.ArgumentParser(description='takes directory of Cv bedgraph files, returns difference peaks within monocistronic CDS')
	parser.add_argument('GFF3',type=str,help='GFF3 annotation file')
	parser.add_argument('monocistronic',type=str,help='list of monocistronic genes')
	parser.add_argument('strand',type=str,help='strand of interest <+/->')

	args = parser.parse_args()

	mono = []
	with open(args.monocistronic,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			mono.append(line.strip())

	annotation = GFF_build(args.GFF3,mono)

	for fyle in sorted(glob('*.bedgraph')):
		if 'BG1_' in fyle:
			BG1 = read_bedgraph(fyle)
		elif 'BG546' in fyle:
			BG546 = read_bedgraph(fyle)
		elif 'BG838' in fyle:
			BG838 = read_bedgraph(fyle)
		elif 'BG1030' in fyle:
			BG1030 = read_bedgraph(fyle)

	BG1_546 = compare(BG1,BG546)
	BG1_838 = compare(BG1,BG838)
	BG1_1030 = compare(BG1,BG1030)
	BG546_838 = compare(BG546,BG838)
	BG546_1030 = compare(BG546,BG1030)
	BG838_1030 = compare(BG838,BG1030)

	penultimate = {
	'BG1_BG546' : BG1_546,
	'BG1_BG838' : BG1_838,
	'BG1_BG1030' : BG1_1030,
	'BG546_BG838' : BG546_838,
	'BG546_BG1030' : BG546_1030,
	'BG838_BG1030' : BG838_1030
	}

	final = {}
	for key,value in penultimate.items():
		thresh = CDS_thresh(value[1],annotation,args.strand)
		final[key] = [value[0],thresh]

	for key,value in final.items():
		with open(key+'.bedgraph','w') as outp:
			for count, subvalue in value[0]:
				if count in value[1]:
					outp.write('NC_000964.3\t'+str(count)+'\t'+str(count+1)+'\t'+str(subvalue)+'\n')

if __name__ == '__main__':
	main()






