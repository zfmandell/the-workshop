#!/usr/bin/env python3

import argparse
import csv

def term_reader(csv):
	return_dict = {}
	with open(csv,'r') as inp:
		firstline = inp.readline()
		for line in inp: 
			coord = int(line.strip().split(',')[0])
			strand = line.strip().split(',')[1]
			dist = int(line.strip().split(',')[2])
			return_dict[coord] = [strand,dist]
	
	return return_dict

def table_reader(table):
	return_dict = {}
	with open(table,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			coord = int(line.strip().split(',')[0])
			contents = line.strip().split(',')[1:]
			return_dict[coord] = contents
			
	return return_dict,firstline

def genome_yield(fasta_name):
	seq = ''
	with open(fasta_name) as inp:
		for line in inp:
			if line[0] != '>':
				seq = seq+str(line.strip())
	return seq

def reverse_complement(seq):
	#returns reverse complement of a list of nucleotides
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	bases = list(seq)
	bases = reversed([complement.get(base,base) for base in bases])
	bases = ''.join(bases)
	return bases

def complement(seq):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	bases = list(seq)
	bases = [complement.get(base,base) for base in bases]
	bases = ''.join(bases)
	return bases

def table_adjuster(terms,table,genome):
	return_dict,i = {},0
	for coord,value in table.items():
		strand = value[0]
		pos = value[1]
		hairpin = value[3]
		A = value[2]
		U = value[4]
		if strand == '+':
			chunk = genome[coord-100:coord].replace('T','U')
		else:
			chunk = reverse_complement(genome[coord:coord+101]).replace('T','U')

		try:
			if terms[coord][0] == strand:
				dist = terms[coord][1]

				if strand == '+':
					if dist == 0:
						adjusted_coord = coord
						adjusted_pos = pos
						adjusted_U = U
						adjusted_hairpin = hairpin
						adjusted_A = A
					else:
						adjusted_coord = coord+dist
						try:
							adjusted_pos = '+'.join([pos.split("+")[0],str(int(pos.split('+')[1].strip('nt'))+dist)])+'nt'
						except IndexError:
							adjusted_pos = pos
						adjusted_U = U[dist:]+genome[coord:coord+dist].replace('T','U')
						adjusted_hairpin = A[-dist:] + hairpin + U[:dist]
						adjusted_A = chunk[chunk.find(A)-1:chunk.find(A)+8]
				else:
					if dist == 0:
						adjusted_coord = coord
						adjusted_pos = pos
						adjusted_U = U
						adjusted_hairpin = hairpin
						adjusted_A = A
					else:
						adjusted_coord = coord-dist
						try:
							adjusted_pos = '+'.join([pos.split("+")[0],str(int(pos.split('+')[1].strip('nt'))+dist)])+'nt'
						except IndexError:
							adjusted_pos = pos
						adjusted_U = U[dist:]+complement(genome[coord-dist:coord]).replace('T','U')
						adjusted_hairpin = A[-dist:] + hairpin + U[:dist]
						adjusted_A = chunk[chunk.find(A)-1:chunk.find(A)+8]
						
			value[1] = adjusted_pos
			value[2] = adjusted_A
			value[3] = adjusted_hairpin
			value[4] = adjusted_U

			return_dict[coord] = value

		except KeyError:
			pass

	return return_dict
					
def main():
	parser = argparse.ArgumentParser(description='determines 3-OH end idendity by combining Term-seq and RNET-seq information')
	parser.add_argument('ends',type=str,help='list of Term-seq called 3-OH ends, adjusted using RNET-seq info: <.csv>, no default value, must be inputted')
	parser.add_argument('Total_Table',type=str,help='Table S2 from NusG Termination Paper: <.csv>, no default value, must be inputted')
	parser.add_argument('genome',type=str,help='genome to be used to fix Table S2, no default value, must be inputted')

	parser.add_argument('output',type=str,help='name of output file')
   
	args = parser.parse_args()

	terms = term_reader(args.ends)
	genome = genome_yield(args.genome)
	table,firstline = table_reader(args.Total_Table)
	adjusted_table = table_adjuster(terms,table,genome)

	with open(args.output,'w') as outp:
		outp.write(firstline)
		for key,value in adjusted_table.items():
			outp.write(str(key)+','+",".join(list(map(str,value)))+'\n')
	  
if __name__ == '__main__':
	main()