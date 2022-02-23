import sys

A = sys.argv[1]
B = sys.argv[2]

A_terms = []
with open(A,'r') as inp:
	firstline = inp.readline()
	for line in inp:
