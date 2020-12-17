#!/usr/bin/env python

import sys


with open(sys.argv[1],'r') as inp:
    firstline,master_dict = inp.readline(),{}
    for line in inp:
        if len(line.strip().split('\t')) != 1:
            master_dict[int(line.strip().split('\t')[1])] = [line.strip().split('\t')[8],line.strip().split('\t')[9]]
        else:
            pass

with open(sys.argv[2],'r') as inp:
    firstline,probe_dict = inp.readline(),{}
    for line in inp:
        if len(line.strip().split(',')) != 1:
            if len(line.strip().split(',')[0].split("_")) == 2:
                probe_dict[int(line.strip().split(',')[0].split("_")[1])] = [line.strip().split(',')[1],line.strip().split(',')[2]]
            else:
                pass
        else:
            pass

return_dict = {}
for key,value in probe_dict.iteritems():
    if 'NA' not in value:
        return_dict[key] = [str(value[0]),str(value[1]),str(abs(float(master_dict[key][0]))),str(abs(float(master_dict[key][1])))]
    else:
        return_dict[key] = [str(value[0]),'NA',str(abs(float(master_dict[key][0]))),str(abs(float(master_dict[key][1])))]

with open(sys.argv[3],'w') as outp:
    outp.write('coord,deltaG,pos,count,score\n')
    for key,value in return_dict.iteritems():
        print value
        outp.write(str(key)+","+value[0]+','+value[1]+','+value[2]+','+value[3]+'\n')
