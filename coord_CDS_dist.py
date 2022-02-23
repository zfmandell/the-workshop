#!/usr/bin/env python3

import argparse

def read_ends(fyle):
	return_list = []
	with open(fyle,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			return_list.append(int(line.strip().split(',')[0]))

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

def dist_calc(peak,start,end,strand):
	if strand == '+':
		return (peak-start)/(end-start+1)
	else:
		return (end-peak)/(end-start+1)

def main():
	parser = argparse.ArgumentParser(description='takes directory of Cv bedgraph files, returns difference peaks within monocistronic CDS')
	parser.add_argument('GFF3',type=str,help='GFF3 annotation file')
	parser.add_argument('monocistronic',type=str,help='list of monocistronic genes')
	parser.add_argument('ends',type=str,help='list of ends <.csvs>')

	args = parser.parse_args()

	mono = []
	with open(args.monocistronic,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			mono.append(line.strip())

	annotation = GFF_build(args.GFF3,mono)
	ends = read_ends(args.ends)

	final = {}
	for item in ends:
		for key,value in annotation.items():
			if value[1] <= item <= value[2]:
				final[item] = dist_calc(item,value[1],value[2],value[0])

	for key,value in final.items():
		print(key,round(value,2))

if __name__ == '__main__':
	main()