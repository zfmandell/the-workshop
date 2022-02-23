#!/usr/bin/env python3

import sys

#removes terms from a bedgraph, makes new bedgraph without terms

terminators = sys.argv[1]
bedgraph = sys.argv[2]

terms = {}
with open(terminators,'r') as inp:
	firstline = inp.readline()
	for line in inp:
		coord = int(line.strip().split(',')[0])
		strand = line.strip().split(',')[1]
		terms[coord] = strand

marked = []
with open(bedgraph,'r') as inp:
	with open('term_free_'+bedgraph,'w') as outp:
		for line in inp:
			coord = int(line.strip().split()[2])
			delta = float(line.strip().split()[3])
			if delta > 0:
				strand = '+'
			else:
				strand = '-'
			for key in terms.keys():
				if key-3 <= coord <= key+3:
					if strand == terms[key]:
						marked.append(coord)
		inp.seek(0)
		for line in inp:
			coord = int(line.strip().split()[2])
			if coord in marked:
				pass
			else:
				outp.write(line)


