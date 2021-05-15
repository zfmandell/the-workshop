#!/usr/bin/env python3

import argparse
from glob import glob
from collections import OrderedDict

def read_bedgraph(bed_fyle):
	#create list of coverage values, position dependent
	return_dict = {}
	with open(bed_fyle,"r") as inp:
		for line in inp:
			coord = int(line.strip().split()[2])
			delta = float(line.strip().split()[3])
			if delta > 0:
				strand = '+'
			else:
				strand = '-'
			return_dict[coord] = [delta,strand]

	return return_dict

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


def CDS_thresh(dyct,GFF3):
	return_dict = {}
	for key,value in dyct.items():
		if queried(key,GFF3,value[1]) == True:
			return_dict[key] = value

	return return_dict

def queried(number,GFF3,strand):
	for key,value in GFF3.items():
		if (strand == value[0]) and (value[1] <= number <= value[2]):
			return True

def main():
	parser = argparse.ArgumentParser(description='takes a Cv files, returns Cv peaks thresholded by monocistronic CDS')
	parser.add_argument('bedgraph',type=str,help='bedgraph file to operate on')
	parser.add_argument('GFF3',type=str,help='GFF3 annotation file')
	parser.add_argument('monocistronic',type=str,help='list of monocistronic genes')

	args = parser.parse_args()

	mono = []
	with open(args.monocistronic,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			mono.append(line.strip())

	annotation = GFF_build(args.GFF3,mono)

	peaks = read_bedgraph(args.bedgraph)
	thresh = CDS_thresh(peaks,annotation)

	with open('mono_CDS_only_'+args.bedgraph,'w') as outp:
		for key,value in OrderedDict(sorted(thresh.items())).items():
			outp.write('NC_000964.3\t'+str(key-1)+'\t'+str(key)+'\t'+str(value[0])+'\n')

if __name__ == '__main__':
	main()
