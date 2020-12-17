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
        base_seq = genome['fasta'][end-200:end+200]
        if seq_find(base_seq,hairpin) == -1:
            print 'no'

        else:
            end = seq_find(base_seq,hairpin)
            start = seq_find(base_seq,hairpin)-12
            return base_seq[start:end].replace('T','U')
    else:
        base_seq = reverse_complement(genome['fasta'][end-200:end+200])
        if seq_find(base_seq,hairpin) == -1:
            print 'no'
        else:
            end = seq_find(base_seq,hairpin)
            start = seq_find(base_seq,hairpin)-12
            return base_seq[start:end].replace('T','U')

def downstream_find(genome,hairpin,end,strand):
    if strand == '+':
        base_seq = genome['fasta'][end-200:end+200]
        if seq_find(base_seq,hairpin) == -1:
            print 'no'
        else:
            end = seq_find(base_seq,hairpin)+len(hairpin)+12
            start = seq_find(base_seq,hairpin)+len(hairpin)
            return base_seq[start:end]
    else:
        base_seq = reverse_complement(genome['fasta'][end-200:end+200])
        if seq_find(base_seq,hairpin) == -1:
            print 'no'
        else:
            end = seq_find(base_seq,hairpin)+len(hairpin)+12
            start = seq_find(base_seq,hairpin)+len(hairpin)
            return base_seq[start:end]

genome = read_fasta(sys.argv[1])
hairpin,up,down = [],[],[]
with open(sys.argv[2],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        hairpin.append([line.strip().split(",")[4],line.strip().split(",")[0],line.strip().split(",")[1]])

for item in hairpin:
    up.append(upstream_find(genome,item[0],int(item[1]),item[2]))
    down.append(downstream_find(genome,item[0],int(item[1]),item[2]))

perc_up,perc_down = {},{}
i = 0
while i < 12:
    temp_up = Counter([x[i] for x in up])
    temp_down = Counter([x[i] for x in down])

    to_add_up,to_add_down = {},{}

    for key,value in temp_up.iteritems():
        to_add_up[key] = (value/len(hairpin))*100

    for key,value in temp_down.iteritems():
        to_add_down[key] = (value/len(hairpin))*100

    perc_up[i] = to_add_up
    perc_down[i] = to_add_down
    i+=1

with open('up_'+sys.argv[2],'w') as outp:
    outp.write('pos,nucleotide,perc\n')
    i = 0
    while i < 12:
        outp.write(str(i-12)+',A,')
        if 'A' in perc_up[i].keys():
            outp.write(str(perc_up[i]['A'])+"\n")
        else:
            outp.write('0\n')
        outp.write(str(i-12)+',U,')
        if 'U' in perc_up[i].keys():
            outp.write(str(perc_up[i]['U'])+"\n")
        else:
            outp.write('0\n')
        outp.write(str(i-12)+',C,')
        if 'C' in perc_up[i].keys():
            outp.write(str(perc_up[i]['C'])+"\n")
        else:
            outp.write('0\n')
        outp.write(str(i-12)+',G,')
        if 'G' in perc_up[i].keys():
            outp.write(str(perc_up[i]['G'])+"\n")
        else:
            outp.write('0\n')
        i += 1

with open('down_'+sys.argv[2],'w') as outp:
    outp.write('pos,nucleotide,perc\n')
    i = 0
    while i < 12:
        outp.write(str(i+1)+',A,')
        if 'A' in perc_down[i].keys():
            outp.write(str(perc_down[i]['A'])+"\n")
        else:
            outp.write('0\n')
        outp.write(str(i+1)+',T,')
        if 'T' in perc_down[i].keys():
            outp.write(str(perc_down[i]['T'])+"\n")
        else:
            outp.write('0\n')
        outp.write(str(i+1)+',C,')
        if 'C' in perc_down[i].keys():
            outp.write(str(perc_down[i]['C'])+"\n")
        else:
            outp.write('0\n')
        outp.write(str(i+1)+',G,')
        if 'G' in perc_down[i].keys():
            outp.write(str(perc_down[i]['G'])+"\n")
        else:
            outp.write('0\n')
        i += 1
