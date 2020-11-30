import sys

with open(sys.argv[1],'r') as inp:
    samples = {}
    for line in inp:
        contents = line.strip().split("\t")
        sample_contents = contents[0].split(".")
        sample = sample_contents[0]
        name = str(sample)+'.txt'

        with open(name,"a+") as outp:
            outp.write(line)
