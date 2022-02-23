#!/usr/bin/env python3

import argparse

def read_bedgraph(bed_fyle):
	#create list of coverage values, position dependent
	return_dict = {}
	with open(bed_fyle,"r") as inp:
		for line in inp:
			coord = int(line.strip().split()[2])
			delta = float(line.strip().split()[3])
			if delta > 0:
				strand = '+'
			else:
				strand = '-'
			return_dict[coord] = [delta,strand]

	return return_dict

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

def seq_strip(sequence,pos,length,strand):
	if strand == '+':
		return sequence[pos-(length):pos].replace('T','U')
	else:
		return reverse_complement(sequence[pos-1:pos+length-1]).replace('T','U')

def main():
	parser = argparse.ArgumentParser(description='takes a Cv file, returns  set of sequences based on peaks')
	parser.add_argument('bedgraph',type=str,help='bedgraph file to operate on')
	parser.add_argument('genome',type=str,help='genomic fasta')
	parser.add_argument('length',type=int,help='length of return sequences')
	parser.add_argument('output',type=str,help='name of output file')

	args = parser.parse_args()

	sequence = genome_yield(args.genome)
	peaks = read_bedgraph(args.bedgraph)

	final = {}
	for key,value in peaks.items():
		if value[1] == '+':
			final['>'+str(key)+','+str(value[0])+','+str(value[1])] = seq_strip(sequence,key,args.length,'+')
		else:
			final['>'+str(key)+','+str(value[0])+','+str(value[1])] = seq_strip(sequence,key,args.length,'-')

	with open(args.output,'w') as outp:
		for key,value in final.items():
			outp.write(key+'\n'+value+"\n")

if __name__ == '__main__':
	main()