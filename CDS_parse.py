#!/usr/bin/env python
import sys
from Bio import SeqIO

for record in SeqIO.parse(sys.argv[1], "fasta"):

    contents = str(record.id).split("|")

    new_contents = []

    for item in contents:
        if ":" in str(item):
            new_contents.append(item.split(":"))
        else:
            new_contents.append(item)

    if str(new_contents[1][0]) == 'CDS':

        sequence = str(record.seq)[0:int(new_contents[1][1])]

        fasta_output = [sequence[i:i+80] for i in range(0, len(sequence), 80)]

        with open(sys.argv[2],'a') as fasta:
            fasta.write(">"+str(record.id))
            fasta.write("\n")
            for item in fasta_output:
                fasta.write(item)
                fasta.write("\n")


    elif str(new_contents[2][0]) == 'CDS' and len(contents) < 4 :

        sequence = str(record.seq)[int(new_contents[1][1]):]

        fasta_output = [sequence[i:i+80] for i in range(0, len(sequence), 80)]

        with open(sys.argv[2],'a') as fasta:
            fasta.write(">"+str(record.id))
            fasta.write("\n")
            for item in fasta_output:
                fasta.write(item)
                fasta.write("\n")

    elif str(new_contents[2][0]) == 'CDS' and len(contents) == 4 :

        sequence = str(record.seq)[int(new_contents[1][1]):(int(new_contents[1][1])+int(new_contents[2][1]))]

        fasta_output = [sequence[i:i+80] for i in range(0, len(sequence), 80)]

        with open(sys.argv[2],'a') as fasta:
            fasta.write(">"+str(record.id))
            fasta.write("\n")
            for item in fasta_output:
                fasta.write(item)
                fasta.write("\n")
