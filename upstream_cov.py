#!/usr/bin/env python3

import sys
from glob import glob
import numpy as np

def read_terms(fyle):
	return_list = []
	with open(fyle,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			coord = line.strip().split(',')[0]
			strand = line.strip().split(',')[1]

			return_list.append(coord+','+strand)
	
	return return_list

def read_coverage(cov_fyle):
	#create list of coverage values, position dependent
	return_list = []
	with open(cov_fyle,"r") as inp:
		for line in inp:
			return_list.append(float(line.strip().split('\t')[2]))
	return return_list

def upstream_calc(coord,cov,strand):
	if strand == '+':
		return np.median(cov[coord-10:coord])
	
	else:
		return np.median(cov[coord+1:coord+11])

def read_cov(fwd,rev,terms):
	WT,A,G,R,NGN,PNP,AG,GR,AR,NGNR,AGR = {},{},{},{},{},{},{},{},{},{},{}
	
	for fyle in sorted(glob(fwd+'/*.cov')):
		print(fyle)
		if fyle.split('/')[-1].split('_')[0] == 'WT':
			temp = read_coverage(fyle)
			for item in terms:
				key = int(item.split(',')[0])
				value = item.split(',')[1]
				if value == '+':
					WT[item] = upstream_calc(key,temp,value)
		elif fyle.split('/')[-1].split('_')[0] == 'A':
			temp = read_coverage(fyle)
			for item in terms:
				key = int(item.split(',')[0])
				value = item.split(',')[1]
				if value == '+':
					A[item] = upstream_calc(key,temp,value)
		elif fyle.split('/')[-1].split('_')[0] == 'G':
			temp = read_coverage(fyle)
			for item in terms:
				key = int(item.split(',')[0])
				value = item.split(',')[1]
				if value == '+':
					G[item] = upstream_calc(key,temp,value)
		elif fyle.split('/')[-1].split('_')[0] == 'R':
			temp = read_coverage(fyle)
			for item in terms:
				key = int(item.split(',')[0])
				value = item.split(',')[1]
				if value == '+':
					R[item] = upstream_calc(key,temp,value)
		elif fyle.split('/')[-1].split('_')[0] == 'NGN':
			temp = read_coverage(fyle)
			for item in terms:
				key = int(item.split(',')[0])
				value = item.split(',')[1]
				if value == '+':
					NGN[item] = upstream_calc(key,temp,value)
		elif fyle.split('/')[-1].split('_')[0] == 'PNP':
			temp = read_coverage(fyle)
			for item in terms:
				key = int(item.split(',')[0])
				value = item.split(',')[1]
				if value == '+':
					PNP[item] = upstream_calc(key,temp,value)
		elif fyle.split('/')[-1].split('_')[0] == 'AG':
			temp = read_coverage(fyle)
			for item in terms:
				key = int(item.split(',')[0])
				value = item.split(',')[1]
				if value == '+':
					AG[item] = upstream_calc(key,temp,value)
		elif fyle.split('/')[-1].split('_')[0] == 'GR':
			temp = read_coverage(fyle)
			for item in terms:
				key = int(item.split(',')[0])
				value = item.split(',')[1]
				if value == '+':
					GR[item] = upstream_calc(key,temp,value)
		elif fyle.split('/')[-1].split('_')[0] == 'AR':
			temp = read_coverage(fyle)
			for item in terms:
				key = int(item.split(',')[0])
				value = item.split(',')[1]
				if value == '+':
					AR[item] = upstream_calc(key,temp,value)
		elif fyle.split('/')[-1].split('_')[0] == 'NGNR':
			temp = read_coverage(fyle)
			for item in terms:
				key = int(item.split(',')[0])
				value = item.split(',')[1]
				if value == '+':
					NGNR[item] = upstream_calc(key,temp,value)
		elif fyle.split('/')[-1].split('_')[0] == 'AGR':
			temp = read_coverage(fyle)
			for item in terms:
				key = int(item.split(',')[0])
				value = item.split(',')[1]
				if value == '+':
					AGR[item] = upstream_calc(key,temp,value)

	for fyle in sorted(glob(rev+'/*.cov')):
		print(fyle)
		if fyle.split('/')[-1].split('_')[0] == 'WT':
			temp = read_coverage(fyle)
			for item in terms:
				key = int(item.split(',')[0])
				value = item.split(',')[1]
				if value == '-':
					WT[item] = upstream_calc(key,temp,value)
		elif fyle.split('/')[-1].split('_')[0] == 'A':
			temp = read_coverage(fyle)
			for item in terms:
				key = int(item.split(',')[0])
				value = item.split(',')[1]
				if value == '-':
					A[item] = upstream_calc(key,temp,value)
		elif fyle.split('/')[-1].split('_')[0] == 'G':
			temp = read_coverage(fyle)
			for item in terms:
				key = int(item.split(',')[0])
				value = item.split(',')[1]
				if value == '-':
					G[item] = upstream_calc(key,temp,value)
		elif fyle.split('/')[-1].split('_')[0] == 'R':
			temp = read_coverage(fyle)
			for item in terms:
				key = int(item.split(',')[0])
				value = item.split(',')[1]
				if value == '-':
					R[item] = upstream_calc(key,temp,value)
		elif fyle.split('/')[-1].split('_')[0] == 'NGN':
			temp = read_coverage(fyle)
			for item in terms:
				key = int(item.split(',')[0])
				value = item.split(',')[1]
				if value == '-':
					NGN[item] = upstream_calc(key,temp,value)
		elif fyle.split('/')[-1].split('_')[0] == 'PNP':
			temp = read_coverage(fyle)
			for item in terms:
				key = int(item.split(',')[0])
				value = item.split(',')[1]
				if value == '-':
					PNP[item] = upstream_calc(key,temp,value)
		elif fyle.split('/')[-1].split('_')[0] == 'AG':
			temp = read_coverage(fyle)
			for item in terms:
				key = int(item.split(',')[0])
				value = item.split(',')[1]
				if value == '-':
					AG[item] = upstream_calc(key,temp,value)
		elif fyle.split('/')[-1].split('_')[0] == 'GR':
			temp = read_coverage(fyle)
			for item in terms:
				key = int(item.split(',')[0])
				value = item.split(',')[1]
				if value == '-':
					GR[item] = upstream_calc(key,temp,value)
		elif fyle.split('/')[-1].split('_')[0] == 'AR':
			temp = read_coverage(fyle)
			for item in terms:
				key = int(item.split(',')[0])
				value = item.split(',')[1]
				if value == '-':
					AR[item] = upstream_calc(key,temp,value)
		elif fyle.split('/')[-1].split('_')[0] == 'NGNR':
			temp = read_coverage(fyle)
			for item in terms:
				key = int(item.split(',')[0])
				value = item.split(',')[1]
				if value == '-':
					NGNR[item] = upstream_calc(key,temp,value)
		elif fyle.split('/')[-1].split('_')[0] == 'AGR':
			temp = read_coverage(fyle)
			for item in terms:
				key = int(item.split(',')[0])
				value = item.split(',')[1]
				if value == '-':
					AGR[item] = upstream_calc(key,temp,value)

	return WT,A,G,R,NGN,PNP,AG,GR,AR,NGNR,AGR

terms = read_terms(sys.argv[3])
WT,A,G,R,NGN,PNP,AG,GR,AR,NGNR,AGR = read_cov(sys.argv[1],sys.argv[2],terms)

printed = {}
for item in terms:
	printed[item] = [item,WT[item],A[item],G[item],R[item],NGN[item],PNP[item],AG[item],GR[item],AR[item],NGNR[item],AGR[item]]

with open(sys.argv[4],'w') as outp:
	outp.write('coord,strand,WT,A,G,R,NGN,PNP,AG,GR,AR,NGNR,AGR\n')
	for item in printed.values():
		outp.write(','.join(list(map(str,item)))+'\n')




