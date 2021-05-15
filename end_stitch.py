#!/usr/bin/env python3

import argparse
from glob import glob
import numpy as np
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
	return_dict = {}
	with open(bed_fyle,"r") as inp:
		for line in inp:
			coord = int(line.strip().split()[2])
			delta = float(line.strip().split()[3])
			if delta > 0:
				strand = '+'
			else:
				strand = '-'
			return_dict[coord] = [strand,delta]
	return return_dict

def read_delta(bed_fyle):
	return_list = []
	with open(bed_fyle,"r") as inp:
		for line in inp:
			return_list.append(float(line.strip().split()[3]))
	return return_list

def end_annotater(end,strand,GFF3):
	for key,value in GFF3.items():
		if (value[0] == strand) and (value[1] <= end<= value[2]):
			return [key,dist_calc(end,value[1],value[2],strand)]

def dist_calc(peak,start,end,strand):
	if strand == '+':
		return round((peak-start)/(end-start+1)*100,2)
	else:
		return round((end-peak)/(end-start+1)*100,2)
	
def reference_build(base,query):
	found = {}
	for item in query.keys():
		for subitem in base.keys():
			if subitem-5 <= item <= subitem+5:
				found[subitem] = item

	return found

def build_total(WT_dyct,dyct1,dyct1r,dyct2,dyct2r,dyct3,dyct3r):
	return_dict = {}
	
	for item in WT_dyct.keys():
		return_dict[item] = WT_dyct[item][0]
	for item in dyct1.keys():
		if item not in dyct1r.values():
			return_dict[item] = dyct1[item][0]
	for item in dyct2.keys():
		if item not in dyct2r.values():
			return_dict[item] = dyct2[item][0]
	for item in dyct3.keys():
		if item not in dyct3r.values():
			return_dict[item] = dyct3[item][0]

	return return_dict

def end_stitch(deltas,dyct1r,dyct2r,dyct3r,total):
	return_dict = {
	'BG1' : {},
	'BG546' : {},
	'BG838' : {},
	'BG1030' : {}
	}

	for fyle in sorted(glob(deltas+'/*')):
	
		if ('BG1_' in fyle) and ('fwd' in fyle):
			temp = read_delta(fyle)
			for key,value in total.items():
				if value == '+':
					return_dict['BG1'][key] = temp[key-1]

		if ('BG1_' in fyle) and ('rev' in fyle):
			temp = read_delta(fyle)
			for key,value in total.items():
				if value == '-':
					return_dict['BG1'][key] = temp[key-1]

		if ('BG546' in fyle) and ('fwd' in fyle):
			temp = read_delta(fyle)
			for key,value in total.items():
				if value == '+':
					if key in dyct1r.keys():
						return_dict['BG546'][key] = temp[dyct1r[key]-1]
					else:
						return_dict['BG546'][key] = temp[key-1]

		if ('BG546' in fyle) and ('rev' in fyle):
			temp = read_delta(fyle)
			for key,value in total.items():
				if value == '-':
					if key in dyct1r.keys():
						return_dict['BG546'][key] = temp[dyct1r[key]-1]
					else:
						return_dict['BG546'][key] = temp[key-1]

		if ('BG838' in fyle) and ('fwd' in fyle):
			temp = read_delta(fyle)
			for key,value in total.items():
				if value == '+':
					if key in dyct2r.keys():
						return_dict['BG838'][key] = temp[dyct2r[key]-1]
					else:
						return_dict['BG838'][key] = temp[key-1]

		if ('BG838' in fyle) and ('rev' in fyle):
			temp = read_delta(fyle)
			for key,value in total.items():
				if value == '-':
					if key in dyct2r.keys():
						return_dict['BG838'][key] = temp[dyct2r[key]-1]
					else:
						return_dict['BG838'][key] = temp[key-1]

		if ('BG1030' in fyle) and ('fwd' in fyle):
			temp = read_delta(fyle)
			for key,value in total.items():
				if value == '+':
					if key in dyct2r.keys():
						return_dict['BG1030'][key] = temp[dyct2r[key]-1]
					else:
						return_dict['BG1030'][key] = temp[key-1]

		if ('BG1030' in fyle) and ('rev' in fyle):
			temp = read_delta(fyle)
			for key,value in total.items():
				if value == '-':
					if key in dyct2r.keys():
						return_dict['BG1030'][key] = temp[dyct2r[key]-1]
					else:
						return_dict['BG1030'][key] = temp[key-1]

	return return_dict

def build_final(total,stitched,annotation):
	return_dict = {}
	for key,value in total.items():
		temp = []
		annotated = end_annotater(key,value,annotation)

		if stitched['BG1'][key] == 0.0:
			stitched['BG1'][key] = 0.01
		
		if stitched['BG546'][key] == 0.0:
			stitched['BG546'][key] = 0.01
		
		if stitched['BG838'][key] == 0.0:
			stitched['BG838'][key] = 0.01
		
		if stitched['BG1030'][key] == 0.0:
			stitched['BG1030'][key] = 0.01

		temp = [value,annotated[0],annotated[1],stitched['BG1'][key],stitched['BG546'][key],stitched['BG838'][key],stitched['BG1030'][key]]
		temp.append(round(np.log2(stitched['BG546'][key]/stitched['BG1'][key]),2))
		temp.append(round(np.log2(stitched['BG838'][key]/stitched['BG1'][key]),2))
		temp.append(round(np.log2(stitched['BG1030'][key]/stitched['BG1'][key]),2))
		temp.append(round(np.std([stitched['BG1'][key],stitched['BG546'][key],stitched['BG838'][key],stitched['BG1030'][key]]),2))

		return_dict[key] = temp

	return return_dict

def main():
	parser = argparse.ArgumentParser(description='takes directory of Cv peak files, delta values at all peaks of all strains')
	parser.add_argument('GFF3',type=str,help='GFF3 annotation file')
	parser.add_argument('monocistronic',type=str,help='list of monocistronic genes')
	parser.add_argument('deltas',type=str,help='directory containing delta files')
	parser.add_argument('output',type=str,help='name out output file')


	args = parser.parse_args()

	mono = []
	with open(args.monocistronic,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			mono.append(line.strip())

	annotation = GFF_build(args.GFF3,mono)

	for fyle in sorted(glob('thresh*.bedgraph')):
		if 'BG1_' in fyle:
			BG1_ends = read_bedgraph(fyle)
		elif 'BG546' in fyle:
			BG546_ends = read_bedgraph(fyle)
		elif 'BG838' in fyle:
			BG838_ends = read_bedgraph(fyle)
		elif 'BG1030' in fyle:
			BG1030_ends = read_bedgraph(fyle)

	BG546_ref = reference_build(BG1_ends,BG546_ends)
	BG838_ref = reference_build(BG1_ends,BG838_ends)
	BG1030_ref = reference_build(BG1_ends,BG1030_ends)

	total = build_total(BG1_ends,BG546_ends,BG546_ref,BG838_ends,BG838_ref,BG1030_ends,BG1030_ref)

	stitched = end_stitch(args.deltas,BG546_ref,BG838_ref,BG1030_ref,total)

	final = build_final(total,stitched,annotation)

	with open(args.output,'w') as outp:
		outp.write('coord,strand,gene,dist,BG1,BG546,BG838,BG1030,log2FC(BG546/BG1),log2FC(BG838/BG1),log2FC(BG1030/BG1),std dev\n')
		for key,value in OrderedDict(sorted(final.items())).items():
			outp.write(str(key)+','+','.join(list(map(str,value)))+'\n')

if __name__ == '__main__':
	main()