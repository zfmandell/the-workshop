import sys
from Bio import SeqIO
import re

fastaFile = sys.argv[1]
outpfile = sys.argv[2]

for record in SeqIO.parse(fastaFile, "fasta"):
    with open(outpfile,"a+") as outp:
        sequence = str(record.id)

        contents = sequence.split("|")

        new_contents = []

        for item in contents:
            if ":" in str(item):
                new_contents.append(item.split(":"))
            else:
                new_contents.append(item)

        if  'ins' not in new_contents[0]:
            outp.write(">"+str(record.id))
            outp.write("\n")
            outp.write(str(record.seq))
            outp.write("\n")
        else:
            print 'found'
