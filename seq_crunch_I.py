from __future__ import division
import sys
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def reverse_complement(seq):
    alt_map = {'ins':'0'}
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    for k,v in alt_map.iteritems():
        seq = seq.replace(k,v)
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.iteritems():
        bases = bases.replace(v,k)
    return bases

def seq_find(seq,query):
    return seq.replace('T','U').find(query)

def read_fasta(genome_fasta):
    fasta_sequences,fasta_dict =SeqIO.parse(open(genome_fasta),'fasta'),{}
    for fasta in fasta_sequences:
        fasta_dict['fasta'] = str(fasta.seq)
    return fasta_dict

def upstream_find(genome,hairpin,end,strand):
    if strand == '+':
        base_seq = genome['fasta'][end-700:end+700]
        if seq_find(base_seq,hairpin) == -1:
            print 'no'

        else:
            end = seq_find(base_seq,hairpin)
            start = seq_find(base_seq,hairpin)-500
            return base_seq[start:end].replace('T','U')
    else:
        base_seq = reverse_complement(genome['fasta'][end-700:end+700])
        if seq_find(base_seq,hairpin) == -1:
            print 'no'
        else:
            end = seq_find(base_seq,hairpin)
            start = seq_find(base_seq,hairpin)-500
            return base_seq[start:end].replace('T','U')

def downstream_find(genome,hairpin,end,strand):
    if strand == '+':
        base_seq = genome['fasta'][end-200:end+200]
        if seq_find(base_seq,hairpin) == -1:
            print 'no'
        else:
            end = seq_find(base_seq,hairpin)+len(hairpin)+15
            start = seq_find(base_seq,hairpin)+len(hairpin)
            return base_seq[start:end]
    else:
        base_seq = reverse_complement(genome['fasta'][end-200:end+200])
        if seq_find(base_seq,hairpin) == -1:
            print 'no'
        else:
            end = seq_find(base_seq,hairpin)+len(hairpin)+15
            start = seq_find(base_seq,hairpin)+len(hairpin)
            return base_seq[start:end]

genome = read_fasta(sys.argv[1])
hairpin,up,down = [],[],[]
with open(sys.argv[2],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        hairpin.append([line.strip().split(",")[4],line.strip().split(",")[0],line.strip().split(",")[1]])

for item in hairpin:
    up.append([upstream_find(genome,item[0],int(item[1]),item[2]),item[1],item[2]])
    down.append([downstream_find(genome,item[0],int(item[1]),item[2]),item[1],item[2]])

with open('up_'+sys.argv[2]+'.fa','w') as outp:
    for item in up:
        outp.write(">"+item[1]+"_"+item[2]+"\n")
        outp.write(item[0]+"\n")
with open('down_'+sys.argv[2]+'.fa','w') as outp:
    for item in down:
        outp.write(">"+item[1]+"_"+item[2]+"\n")
        outp.write(item[0]+"\n")
