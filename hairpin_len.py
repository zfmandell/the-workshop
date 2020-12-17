#!/usr/bin/env python

import sys

with open(sys.argv[1],'r') as inp:
    with open(sys.argv[2],'w') as outp:
        firstline = inp.readline()
        outp.write('len\n')
        for line in inp:
            outp.write(str(len(line.strip().split(",")[10]))+"\n")
