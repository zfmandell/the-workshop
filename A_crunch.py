#!/usr/bin/env python

from __future__ import division
import glob
import csv

def A_crunch(fyle):
    count,tab = 0,0
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            if line.strip().split(",")[2][-3:] == 'AAA':
                tab+=1
            count +=1
    return round((tab/count)*100,3)


for fyle in glob.glob('*.csv'):
    print fyle
    print A_crunch(fyle)
