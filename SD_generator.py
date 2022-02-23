import sys

def genome_yield(fasta_name):
    seq = ''
    with open(fasta_name) as inp:
        for line in inp:
            if line[0] != '>':
                seq = seq+str(line.strip())
    return seq

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases

pauses = sys.argv[1]
raw = sys.argv[2]
genome = genome_yield(sys.argv[3])

pause_dict = {}
with open(pauses,'r') as inp:
	firstline = inp.readline()
	for line in inp:
		coord = int(line.strip().split(',')[0])
		strand = line.strip().split(',')[1]
		gene = line.strip().split(',')[2]
		strength = float(line.strip().split(',')[3])

		if gene not in pause_dict.keys():
			pause_dict[gene] = [coord,strand,strength]
		else:
			if abs(pause_dict[gene][2]) > abs(strength):
				pass 
			else:
				pause_dict[gene] = [coord,strand,strength]

total = [x[0] for x in pause_dict.values()]

with open(raw,'r') as inp:
	firstline = inp.readline()
	for line in inp:
		coord = int(line.strip().split(',')[0])
		gene = line.strip().split(',')[2]
		sense = line.strip().split(',')[6]
		start = int(line.strip().split(',')[13])

		if coord in total:
			if sense == 'SENSE':
				pause_dict[gene].append(start)
			else:
				del pause_dict[gene]

for key,value in pause_dict.items():
	if value[1] == '+':
		pause_dict[key].append(genome[value[3]-21:value[3]-5].lower().replace('t','u'))
	else:
		pause_dict[key].append(reverse_complement(genome[value[3]+4:value[3]+20]).lower().replace('t','u'))

for key,value in pause_dict.items():
	print(key,value)

with open('SD.fasta','w') as outp:
	for key,value in pause_dict.items():
		outp.write('>'+key+'\n'+value[4]+'\n')



