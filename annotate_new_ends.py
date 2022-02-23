#!/usr/bin/env python3

import argparse
import sys
from glob import glob

def build_GFF3(fyle):
	genes,rRNA = [],[]
	with open(fyle,'r') as inp:
		for line in inp:
			if len(line.strip().split('\t')) > 1:
				if line.strip().split('\t')[2] == 'gene':
					gene = line.strip().split('\t')[8].split(';')[2].split('=')[1]
					start = int(line.strip().split('\t')[3])
					end = int(line.strip().split('\t')[4])
					strand = line.strip().split('\t')[6]
					genes.append([gene,strand,start,end])
				   
	return genes

def grab_CT(directory):
	return_dict = {}
	for fyle in sorted(glob(directory+'/*.ct')):
		coord = int(fyle.strip().split('/')[3].split('_')[0])
		with open(fyle,'r') as inp:
			firstline = inp.readline()
			dG = firstline.split()[3]
			return_dict[coord] = dG
	
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

def read_terms(fyle):
	terms = {}
	with open(fyle) as inp:
		firstline = inp.readline()
		for line in inp:
			coord = int(line.strip().split(',')[0])
			strand = line.strip().split(',')[1]
			hairpin = line.strip().split(',')[2]
			dG = line.strip().split(',')[3]

			terms[coord] = [strand,hairpin,dG]
	
	return terms

def GFF3_annotate(ends,plus,minus):
	return_dict = {}
	for key,value in ends.items():
		if value[0] == '+':
			for item in plus[:len(plus)-1]:
				if key > item[3] and key < plus[plus.index(item)+1][2]:
					return_dict[key] = item[0]+'+'+str(key-item[3])+'nt'
		else:
			for item in minus[:len(minus)-1]:
				if key < item[2] and key > minus[minus.index(item)+1][3]:
					return_dict[key] = item[0]+'+'+str(item[2]-key)+'nt'

	return return_dict

def genome_annotate(ends,sequence):
	return_dict = {}
	for key,value in ends.items():

		if value[0] == '+':
			hairpin = value[1]
			probe = sequence[key-100:key+15].replace('T','U')
			loc = probe.find(hairpin)
			if loc != -1:
				A = probe[loc-9:loc]
				U = probe[loc+len(hairpin):loc+len(hairpin)+9]
				return_dict[key] = [A,U]
			else:
				print(key)

		else:
			hairpin = value[1]
			probe = reverse_complement(sequence[key-15:key+100]).replace('T','U')
			loc = probe.find(hairpin)
			if loc != -1:
				A = probe[loc-9:loc]
				U = probe[loc+len(hairpin):loc+len(hairpin)+9]
				return_dict[key] = [A,U]
			else:
				print(key)

	return return_dict
			
def main():
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('new_terms',type=str,help='')
	parser.add_argument('CT',type=str,help='')
	parser.add_argument('GFF3',type=str,help='')
	parser.add_argument('fasta',type=str,help='')
	parser.add_argument('output',type=str,help='')

	args = parser.parse_args()

	new_terms = read_terms(args.new_terms)
	#dG = grab_CT(args.CT)
	GFF3 = build_GFF3(args.GFF3)
	genes_plus = [x for x in GFF3 if x[1] == '+']
	genes_minus = [x for x in GFF3 if x[1] == '-'][::-1]
	sequence = genome_yield(args.fasta)
	gene_annotation = GFF3_annotate(new_terms,genes_plus,genes_minus)
	seq_annotation = genome_annotate(new_terms,sequence)
	
	with open(args.output,'w') as outp:
		outp.write('POT,Strand,Relative Position,upstream sequence,predicted terminator hairpin sequence,downstream sequence,dG-Hairpin (kcal/mol)\n')
		for key,value in new_terms.items():
			outp.write(str(key)+','+value[0]+','+gene_annotation[key]+','+seq_annotation[key][0]+','+value[1]+','+seq_annotation[key][1]+','+value[2]+'\n')
	
if __name__ == '__main__':
	main()