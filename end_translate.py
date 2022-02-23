#!/usr/bin/env python3

import sys

base = []
with open(sys.argv[1],'r') as inp:
	for line in inp:
		base.append(int(line.strip().split()[2]))

query,passed = [],[]
with open(sys.argv[2],'r') as inp:
	for line in inp:
		coord = int(line.strip().split()[2])
		for item in base:
			if item-5 <= coord <= item+5:
				query.append(item)
				passed.append(coord)
	inp.seek(0)
	for line in inp:
		coord = int(line.strip().split()[2])
		if coord in passed:
			pass
		else:
			query.append(coord)

for item in query:
	print(item)



