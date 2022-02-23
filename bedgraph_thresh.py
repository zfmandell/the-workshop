#!/usr/bin/env python3

import sys

with open(sys.argv[1],'r') as inp:
	with open('thresh_'+sys.argv[1],'w') as outp:
		for line in inp:
			delta = float(line.strip().split()[3])

			if abs(delta) >= int(sys.argv[2]):
				outp.write(line)

