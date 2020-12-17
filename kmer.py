from __future__ import division
import subprocess
import sys
import operator
from Bio import SeqIO

final = {}

outp_file = sys.argv[1]
fastaFile = sys.argv[2]

for record in SeqIO.parse(fastaFile, "fasta"):
    sequence = str(record.seq)

    length = len(sequence)

    i = 1
    a = 0
    b = 9

    if length >= 9:
        while i <= length-8:
            seq_parse = sequence[a:b]

            if seq_parse in final.keys():
                final[seq_parse] += 1
            else:
                final[seq_parse] = 1

            i += 1
            a += 1
            b += 1

total = 0
averages = {}

for value in final.values():
    total += value

for value in final.items():
    averages[value[0]] = value[1]/total


sorted_kmers = sorted(averages.items(), key=operator.itemgetter(1),reverse=True)

with open(outp_file,"a+") as outp:
    outp.write("kmer,unique_occurance,normalized_occurence")
    outp.write("\n")
    for item in sorted_kmers:
        outp.write(str(item[0]))
        outp.write(",")
        outp.write(str(item[1]*total))
        outp.write(",")
        outp.write(str(item[1]))
        outp.write("\n")
