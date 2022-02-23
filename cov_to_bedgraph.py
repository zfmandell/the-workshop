#!/usr/bin/env python3

import sys

with open(sys.argv[1],'r') as inp:
	with open(sys.argv[2],'w') as outp:
		for line in inp:
			genome = line.strip().split()[0]
			start = line.strip().split()[1]
			end = str(int(line.strip().split()[1])+1)
			cov = line.strip().split()[2]
			outp.write(genome+'\t'+start+'\t'+end+'\t'+cov+'\n')