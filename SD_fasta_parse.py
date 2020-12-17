#!/usr/bin/env python
import sys
from Bio import SeqIO


for record in SeqIO.parse(sys.argv[1], "fasta"):

    output_name = "SD_"+str(sys.argv[1])

    contents = str(record.id).split("|")

    new_contents = []

    for item in contents:
        if ":" in str(item):
            new_contents.append(item.split(":"))
        else:
            new_contents.append(item)

    if str(new_contents[1][0]) == 'FPUTR' and int(new_contents[1][1]) == 40:

        sequence = str(record.seq)[0:51]

        with open(output_name,'a') as fasta:
            fasta.write(">"+str(record.id))
            fasta.write("\n")
            fasta.write(sequence)
            fasta.write("\n")

    elif str(new_contents[1][0]) == 'FPUTR' and int(new_contents[1][1]) > 40:

        distance = int(new_contents[1][1])-40

        sequence = str(record.seq)[distance:distance+51]

        with open(output_name,'a+') as fasta:
            fasta.write(">"+str(record.id))
            fasta.write("\n")
            fasta.write(sequence)
            fasta.write("\n")
