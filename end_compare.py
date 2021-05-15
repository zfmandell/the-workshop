import sys

ends = sys.argv[1]
bed = sys.argv[2]

total_ends = {}
with open(ends,'r') as inp: 
	firstline = inp.readline()
	for line in inp:
		total_ends[int(line.strip().split(',')[0])] = line.strip().split(',')[1]

print('coord')
with open(bed,'r') as inp:
	for line in inp:
		coord = int(line.strip().split()[2])
		delta = float(line.strip().split()[3])
		
		
		if delta > 0:
			strand = '+'
		else:
			strand = '-'
		for key,value in total_ends.items():
			if (value == strand) and (key-5 <= coord <= key+5):
				print(key)