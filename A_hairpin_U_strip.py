#!/usr/bin/env python3

import sys

with open(sys.argv[1],'r') as inp:
	return_dict = {}
	firstline = inp.readline()
	for line in inp:
		coord = line.strip().split(',')[0]
		strand = line.strip().split(',')[1]
		A = line.strip().split(',')[3]
		hairpin = line.strip().split(',')[4]
		U = line.strip().split(',')[5]

		if sys.argv[2] == 'A':
			return_dict[coord+'_'+strand] = A
		elif sys.argv[2] == 'hairpin':
			return_dict[coord+'_'+strand] = hairpin[-8:]
		elif sys.argv[2] == 'U':
			return_dict[coord+'_'+strand] = U
		elif sys.argv[2] == 'hairpin_U':
			return_dict[coord+'_'+strand] = (hairpin[-6:]+U).replace('U','T')

return_list = []
with open(sys.argv[3],'r') as inp:
	firstline = inp.readline()
	for line in inp:
		coord = line.strip().split(',')[0]
		strand = line.strip().split(',')[1]
		return_list.append(coord+'_'+strand)

final_list = []
with open(sys.argv[4],'r') as inp:
	firstline = inp.readline()
	for line in inp:
		coord = line.strip().split(',')[0]
		strand = line.strip().split(',')[1]
		if coord+'_'+strand in return_list:
			final_list.append(coord+'_'+strand)

with open(sys.argv[5],'w') as outp:
	for item in final_list:
		outp.write('>'+item+'\n'+return_dict[item]+'\n')