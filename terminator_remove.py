#!/usr/bin/env python

import sys

terms = {}
with open(sys.argv[1],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        terms[int(line.strip().split(',')[0])] = line.strip().split(',')[1]

ends = {}
with open(sys.argv[2],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        ends[int(line.strip().split(',')[0])] = line.strip().split(',')[1]

with open(sys.argv[3],'w') as outp:
    outp.write('end,strand\n')
    to_remove = []
    for key,value in ends.iteritems():
        for sub_key,sub_value in terms.iteritems():
            if sub_key-3 <= key <= sub_key+3 and value == sub_value:
                to_remove.append(key)

    for item in to_remove:
        del ends[item]

    for key,value in ends.iteritems():
        outp.write(str(key)+','+str(value)+'\n')
