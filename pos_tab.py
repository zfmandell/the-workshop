#!/usr/bin/env python

from __future__ import division
import sys

with open(sys.argv[1],'r') as inp:
    firstline,ten,eleven,twelve,count = inp.readline(),0,0,0,0
    for line in inp:
        if line.strip().split(',')[2] != 'NA':
            if int(line.strip().split(',')[2]) == 10:
                ten+=1
                count+=1
            elif int(line.strip().split(',')[2]) == 11:
                eleven+=1
                count+=1
            elif int(line.strip().split(',')[2]) == 12:
                twelve+=1
                count+=1
            else:
                count+=1
        else:
            count+=1

print ten/count*100
print eleven/count*100
print twelve/count*100
