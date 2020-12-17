#!/usr/bin/env python

import pandas as pd
import sys

fyle = sys.argv[1]
data = pd.read_excel(fyle) #reading file
IDS = []
with open(sys.argv[2]+'.fasta','w') as outp:
    for index,row in data.iterrows():
        IDS.append(str(row['Genex'])+"-"+str(row['Peak_Location'])+"-"+str(row['Scorex']))
        outp.write(">"+str(row['Genex'])+"-"+str(row['Peak_Location'])+"-"+str(row['Scorex'])+"\n")
        outp.write(str(row['RNA'][10:])+'\n')
with open('IDS_'+sys.argv[2]+'.txt','w') as outp:
    for item in IDS:
        outp.write(item+"\n")
