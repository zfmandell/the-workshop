#!/usr/bin/env python3

import sys

As,hairpins,Us = {},{},{}
with open(sys.argv[1],'r') as inp:
	firstline = inp.readline()
	for line in inp:
		item = line.strip().split(',')[0]
		delta = float(line.strip().split(',')[5])
	
		A = line.strip().split(',')[6]
		hairpin = line.strip().split(',')[7].upper()[-4:]
		U = line.strip().split(',')[8].upper()

		As[item] = A
		hairpins[item] = hairpin
		Us[item] = U 
"""
with open("A_"+sys.argv[2],'w') as outp:
	for key,value in As.items():
		outp.write('>'+key+'\n'+value+'\n')
"""
with open("hairpin_"+sys.argv[2],'w') as outp:
	for key,value in hairpins.items():
		outp.write('>'+key+'\n'+value+'\n')
"""
with open("U_"+sys.argv[2],'w') as outp:
	for key,value in Us.items():
		outp.write('>'+key+'\n'+value+'\n')
"""