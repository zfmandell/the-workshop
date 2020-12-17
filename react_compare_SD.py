import sys
from Bio.Seq import Seq
from Bio import SeqIO

with open(sys.argv[1],"rU") as wild_type:
    with open(sys.argv[2],"rU") as knock_out:

        transcript_wt = wild_type.readline()
        transcript_knockout = knock_out.readline()

        transcript_wt = transcript_wt.strip()

        contents = str(transcript_wt).split("|")

        new_contents = []

        for item in contents:
            if ":" in str(item):
                new_contents.append(item.split(":"))
            else:
                new_contents.append(item)

        with open(sys.argv[3],"rU") as transcriptome:

            for record in SeqIO.parse(transcriptome, "fasta"):

                if str(record.id) == str(transcript_wt):

                    sequence = record.seq[int(sys.argv[4]):int(sys.argv[5])]

                    with open(str(new_contents[0])+"_react_compare.csv","a+") as outp:
                        outp.write("wt,ko,sequence")
                        outp.write("\n")

                        i = 0

                        wild_type_contents = wild_type.readline().split("\t")
                        ko_contents = knock_out.readline().split("\t")

                        while i <= len(sequence)-1:

                            outp.write(wild_type_contents[i])
                            outp.write(",")
                            outp.write(ko_contents[i])
                            outp.write(",")
                            outp.write(sequence[i])
                            outp.write("\n")

                            i += 1
