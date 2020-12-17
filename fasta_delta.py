import sys
from Bio.Seq import Seq
from Bio import SeqIO

contents = []
length = []

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

length = file_len(sys.argv[2])

with open(sys.argv[1], "rU") as handle:
    with open(sys.argv[2],"rU") as inp:

        first_line = inp.readline()
        first_line=first_line.strip()

        for record in SeqIO.parse(handle, "fasta"):

            contents = str(record.id).split("|")

            new_contents = []

            for item in contents:
                if ":" in str(item):
                    new_contents.append(item.split(":"))
                else:
                    new_contents.append(item)


            if str(new_contents[0]) == sys.argv[4]:

                sequence = record.seq[0:length-1]

                i = 0

                with open(sys.argv[3],"a+") as outp:

                    outp.write(first_line+",Seq"+"\n")

                    for line in inp:
                        line = line.strip()
                        outp.write(line+","+sequence[i]+"\n")
                        i += 1
