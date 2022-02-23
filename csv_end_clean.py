import sys

base = {}
with open(sys.argv[1],'r') as inp:
	firstline = inp.readline()
	for line in inp:
		coord = int(line.strip().split(',')[0])
		log2FC = line.strip().split(',')[8:11]
		base[coord] = [log2FC,line]

marked = []
for key in base.keys():
	for subkey in base.keys():
		if key == subkey:
			pass
		else:
			if subkey-5 <= key <= subkey+5:
				if max(base[key][0]) > max(base[subkey][0]):
					marked.append(subkey)
				else:
					marked.append(key)

for item in set(marked):
	del base[item]

with open(sys.argv[2],'w') as outp:
	outp.write(firstline)
	for key,value in base.items():
		outp.write(value[1])

