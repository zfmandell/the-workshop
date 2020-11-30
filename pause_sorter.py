#!/usr/bin/env python

import sys
import pandas as pd

fyle = sys.argv[2]
data = pd.read_excel(fyle) #reading file
downstream,upstream = [],[]
for index,row in data.iterrows():
    downstream.append(int(row['Downstream']))
    upstream.append(int(row['Upstream']))

d,u = 0,0
with open(sys.argv[1],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        if int(line.strip().split(",")[1]) in downstream:
            d+=1
        elif int(line.strip().split(",")[1]) in upstream:
            u+=1

print d,u
