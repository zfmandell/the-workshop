#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import operator

def build_end_list(fyle):
    return_list = []
    with open(fyle,'r') as inp:
        for line in inp:
            return_list.append(int(line.strip('\n')))
    return return_list

def build_WT_dict(fyle):
    return_dict = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            if line.strip().split(",")[1] == 'fwd':
                return_dict[int(line.strip().split(",")[0])] = '+'
            else:
                return_dict[int(line.strip().split(",")[0])] = '-'
    return return_dict

def build_duplex_dict(fyle):
    return_dict = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            return_dict[int(line.strip().split(",")[0])] = line.strip().split(",")[2]
    return return_dict

def supplemental_build(fyle):
    return_dict = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            return_dict[int(line.strip().split(",")[0])] = line.strip().split(",")[1]
    return return_dict

def read_fasta(genome_fasta):
    return str(SeqIO.parse(open(genome_fasta,'fasta').seq))

def annotation_build(fyle):
    return_dict = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            return_dict[int(line.strip().split(",")[0])] = line.strip().split(",")[2]
    return return_dict

def deltaG_build(fyle):
    return_dict = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            return_dict[int(line.strip().split(",")[0])] = line.strip().split(",")[1]
    return return_dict

def translator(lyst,base):
    for item in base:
        if item[0]-3 <= lyst[0] <= item[0]+3 and lyst[1] == item[1]:
            return item[0]

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
        base_seq = genome['fasta'][end-400:end+400]
        if seq_find(base_seq,hairpin) == -1:
            print 'no'

        else:
            end = seq_find(base_seq,hairpin)
            start = seq_find(base_seq,hairpin)-9
            return base_seq[start:end].replace('T','U')
    else:
        base_seq = reverse_complement(genome['fasta'][end-400:end+400])
        if seq_find(base_seq,hairpin) == -1:
            print 'no'
        else:
            end = seq_find(base_seq,hairpin)
            start = seq_find(base_seq,hairpin)-9
            return base_seq[start:end].replace('T','U')

def downstream_find(genome,hairpin,end,strand):
    if strand == '+':
        base_seq = genome['fasta'][end-200:end+200]
        if seq_find(base_seq,hairpin) == -1:
            print 'no'
        else:
            end = seq_find(base_seq,hairpin)+len(hairpin)+9
            start = seq_find(base_seq,hairpin)+len(hairpin)
            return base_seq[start:end].replace('T','U')
    else:
        base_seq = reverse_complement(genome['fasta'][end-200:end+200])
        if seq_find(base_seq,hairpin) == -1:
            print 'no'
        else:
            end = seq_find(base_seq,hairpin)+len(hairpin)+9
            start = seq_find(base_seq,hairpin)+len(hairpin)
            return base_seq[start:end].replace('T','U')

def annotation_sup_build(fyle):
    return_dict = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            return_dict[int(line.strip().split(",")[0])] = line.strip().split(",")[2]
    return return_dict

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('ends',type=str,help='')
    parser.add_argument('WT',type=str,help='')
    parser.add_argument('duplex_old',type=str,help='')
    parser.add_argument('duplex_new',type=str,help='')
    parser.add_argument('supplemental',type=str,help='')
    parser.add_argument('genome',type=str,help='')
    parser.add_argument('annotation',type=str,help='')
    parser.add_argument('deltaG',type=str,help='')
    parser.add_argument('outfile',type=str,help='')
    args = parser.parse_args()

    ends,WT,supplemental = build_end_list(args.ends),build_WT_dict(args.WT),supplemental_build(args.supplemental)
    duplex_new,duplex_old = build_duplex_dict(args.duplex_new),build_duplex_dict(args.duplex_old)
    genome,annotation,deltaG = read_fasta(args.genome),annotation_build(args.annotation),deltaG_build(args.deltaG)
    annotation_sup = annotation_sup_build(args.supplemental)

    ends_translated = {}
    for item in ends:
        if item in duplex_new.keys():
            ends_translated[item] = None
        else:
            ends_translated[item] = translator([item,WT[item]],supplemental.items())

    return_dict = {}
    for key,value in WT.iteritems():
        if (key in ends) and (key in duplex_new.keys()):
            return_dict[key] = [value,annotation[key],upstream_find(genome,duplex_new[key],key,value),duplex_new[key].replace('T','U'),downstream_find(genome,duplex_new[key],key,value),deltaG[key]]

        elif (key in ends) and (key not in duplex_new.keys()):
            return_dict[ends_translated[key]] = [value,annotation_sup[ends_translated[key]],upstream_find(genome,duplex_old[ends_translated[key]],key,value),duplex_old[ends_translated[key]].replace('T','U'),downstream_find(genome,duplex_old[ends_translated[key]],key,value),deltaG[ends_translated[key]]]

    with open(args.outfile,'w') as outp:
        outp.write('POT,strand,relative position,upstream sequence,hairpin sequence,U-tract,G (kcal/mol) of Hairpin\n')
        for item in sorted(return_dict.items(), key=operator.itemgetter(0)):
            outp.write(str(item[0])+','+str(item[1][0])+','+str(item[1][1])+','+str(item[1][2])+','+str(item[1][3])+','+str(item[1][4])+','+str(item[1][5])+'\n')

if __name__ == '__main__':
    main()
