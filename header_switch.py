import sys
from Bio import SeqIO

fasta_1 = {}
fasta_2 = {}
fasta_3 = {}

for record in SeqIO.parse(sys.argv[1], "fasta"):
    sequence = record.seq
    header = record.id

    fasta_1[str(header)] = str(sequence)

for record in SeqIO.parse(sys.argv[2], "fasta"):
    #print record.id
    #sequence = record.seq
    header = record.id

    fasta_2[str(header)] = str(sequence)

for key_1, value_1 in fasta_1.iteritems():


    contents_1 = key_1.split("|")

    transcript_1 = contents_1[0]

    for key_2, value_2 in fasta_2.iteritems():

        contents_2 = key_2.split("::")

        transcript_2 = contents_2[0]

        if str(transcript_1) == str(transcript_2):

            fasta_3[key_2] = value_1



with open("test_fasta.fa","a+") as outp:
    for key, value in fasta_3.iteritems():

        outp.write(">")
        outp.write(key)
        outp.write("\n")
        outp.write(value)
        outp.write("\n")
