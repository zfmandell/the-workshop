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

genome = genome_yield(sys.argv[1])

strand = sys.argv[2]

if strand == '+':
    print(genome[int(sys.argv[3]):int(sys.argv[4])+len(sys.argv[5])])
else:
    print(reverse_complement(genome[int(sys.argv[3]):int(sys.argv[4])+len(sys.argv[5])]))
