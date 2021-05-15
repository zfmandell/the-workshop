#!/usr/bin/env python3

import sys
import glob
from collections import OrderedDict

delta = {}
for fyle in glob.glob('*.ct'):
	with open(fyle,'r') as inp:
		firstline = inp.readline()
		try:
			dG = str(firstline.strip().split()[3])
			coord = fyle.strip().split(',')[0]
			strand = fyle.strip().split(',')[1].split('_')[0]
			
		except IndexError:
			dG = '0'
			coord = fyle.strip().split(',')[0]
			strand = fyle.strip().split(',')[1].split('_')[0]
		
		delta[coord] = [strand,dG]

with open(sys.argv[1],'w') as outp:
	outp.write('POT,strand,deltaG\n')
	for key,value in OrderedDict(sorted(delta.items())).items():
		outp.write(key+','+value[0]+','+value[1]+"\n")
