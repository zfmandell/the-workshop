import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import math

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
            start = seq_find(base_seq,hairpin)-50
            return base_seq[start:end].replace('T','U')
    else:
        base_seq = reverse_complement(genome['fasta'][end-400:end+400])
        if seq_find(base_seq,hairpin) == -1:
            print 'no'
        else:
            end = seq_find(base_seq,hairpin)
            start = seq_find(base_seq,hairpin)-50
            return base_seq[start:end].replace('T','U')

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('ends',type=str,help='')
    parser.add_argument('genome',type=str,help='')
    parser.add_argument('outfile',type=str,help='')
    args = parser.parse_args()

    genome = read_fasta(args.genome)

    with open(args.ends,'r') as inp:
        with open(args.outfile,'w') as outp:
            firstline = inp.readline()
            for line in inp:
                outp.write('>'+line.strip().split(",")[0]+"_"+line.strip().split(",")[1]+"\n")
                hairpin = line.strip().split(",")[4]
                upstream = upstream_find(genome,hairpin,int(line.strip().split(",")[0]),line.strip().split(",")[1])
                compete = hairpin[:int(math.ceil(len(hairpin)*0.75))]
                total = upstream+compete
                outp.write(total+"\n")


if __name__ == '__main__':
    main()
