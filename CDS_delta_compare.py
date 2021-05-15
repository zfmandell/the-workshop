#!/usr/bin/env python3

import argparse
import numpy as np
from glob import glob
from scipy.signal import argrelextrema
from collections import OrderedDict

def read_coverage(cov_fyle):
	#create list of coverage values, position dependent
	return_list = []
	with open(cov_fyle,"r") as inp:
		for line in inp:
			return_list.append(float(line.strip().split('\t')[2]))
	return return_list

def TPM_scaling(CDS,cov):
	RPK = []
	for key,value in CDS.items():
		RPK.append(sum(cov[value[1]:value[2]+1])/(abs(value[1]-value[2])/1000))
	return sum(RPK)/1000000

def TPM_calc(little,big,cov,scaling):
	RPK = sum(cov[little:big+1])/(abs(little-big)/1000)
	return RPK/scaling
	
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

def CDS_crunch(GFF3,fyle,scaling,cov):
	return_dict,coords = {},{}

	with open(fyle,'r') as inp:
		for line in inp:
			delta = float(line.strip().split()[3])  
			if delta != 0:
				coords[int(line.strip().split()[2])] = float(line.strip().split()[3])

	for key,value in GFF3.items():
		TPM = TPM_calc(value[1],value[2],cov,scaling)
		for item in range(value[1],value[2]+1):
			if item in coords.keys():
				if coords[item] != 0 and TPM > 1:
					return_dict[item] = (abs(coords[item]))
					#return_dict[item] = (abs(coords[item]/np.log10(TPM)))
				else:
					return_dict[item] = 0
			else:
				return_dict[item] = 0 

	return return_dict

def CDS_compare(base,query):
	difference,return_list = {},[]
	for key,value in base.items():
		difference[key] = value-query[key]
	
	ordered = OrderedDict(sorted(difference.items()))
	ordered_list = list(ordered)
	arrayed = np.asarray([x for x in ordered.values()])
	maxima = argrelextrema(arrayed, np.greater)[0].tolist()
	minima = argrelextrema(arrayed, np.less)[0].tolist()
	
	for item in maxima:
		return_list.append('NC_000964.3\t'+str(ordered_list[item]-1)+'\t'+str(ordered_list[item])+'\t'+str(difference[ordered_list[item]]))

	return return_list
		
def main():
	parser = argparse.ArgumentParser(description='k-means clustering of monocistronic CDS based on strains BG1,BG546,BG838')
	parser.add_argument('bedgraph',type=str,help='directory of CDS only bedgraph files')
	parser.add_argument('cov',type=str,help='directory of coverage files')
	parser.add_argument('GFF3',type=str,help='GFF3 annotation file')
	parser.add_argument('monocistronic',type=str,help='list of monocistronic genes')
	parser.add_argument('output',type=str,help='name of output std dev file')

	args = parser.parse_args()

	mono = []
	with open(args.monocistronic,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			mono.append(line.strip())

	annotation = GFF_build(args.GFF3,mono)

	for fyle in sorted(glob(args.cov+'/*.cov')):
		if 'BG1_' in fyle:
			BG1_cov = read_coverage(fyle)
		elif 'BG546' in fyle:
			BG546_cov = read_coverage(fyle)
		elif 'BG838' in fyle:
			BG838_cov = read_coverage(fyle)
		elif 'BG1030' in fyle:
			BG1030_cov = read_coverage(fyle)
	
	scaling_BG1 = TPM_scaling(annotation,BG1_cov)
	scaling_BG546 = TPM_scaling(annotation,BG546_cov)
	scaling_BG838 = TPM_scaling(annotation,BG838_cov)
	scaling_BG1030 = TPM_scaling(annotation,BG1030_cov)
	
	for fyle in sorted(glob(args.bedgraph+'/CDS_plus*.bedgraph')):
		if 'BG1.' in fyle:
			CDS_BG1 = CDS_crunch(annotation,fyle,scaling_BG1,BG1_cov)
			BG1_cov = []
		elif 'BG546.' in fyle:
			CDS_BG546 = CDS_crunch(annotation,fyle,scaling_BG546,BG546_cov)
			BG546_cov = []
		elif 'BG838.' in fyle:
			CDS_BG838 = CDS_crunch(annotation,fyle,scaling_BG838,BG838_cov)
			BG838_cov = []
		elif 'BG1030.' in fyle:
			CDS_BG1030 = CDS_crunch(annotation,fyle,scaling_BG1030,BG1030_cov)
			BG1030_cov = []

	one_546 = CDS_compare(CDS_BG1,CDS_BG546)
	one_838 = CDS_compare(CDS_BG1,CDS_BG838)
	one_1030 = CDS_compare(CDS_BG1,CDS_BG1030)
	fivefoursix_838 = CDS_compare(CDS_BG546,CDS_BG838)
	fivefoursix_1030 = CDS_compare(CDS_BG546,CDS_BG1030)
	eightthreeeight_1030 = CDS_compare(CDS_BG838,CDS_BG1030)

	final = {
	'BG1_BG546' : one_546,
	'BG1_BG838' : one_838,
	'BG1_BG1030' : one_1030,
	'BG546_BG838' : fivefoursix_838,
	'BG546_BG1030' : fivefoursix_1030,
	'BG838_BG1030' : eightthreeeight_1030
	}

	for key,value in final.items():
		with open(key+'.bedgraph','w') as outp:
			for subitem in value:
				outp.write(subitem+'\n')
		
		
if __name__ == '__main__':
	main()




