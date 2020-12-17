#!/usr/bin/env python

from __future__ import division
import argparse
from glob import glob
import csv

def deltaT_dict_build(fyle):
    return_dict = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            return_dict[int(line.strip().split(",")[0])] = float(line.strip().split(",")[2])
    return return_dict


def seq_strip(fyle):
    return_dict = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            return_dict[int(line.strip().split(",")[0])] =[str(line.strip().split(",")[2]),str(line.strip().split(",")[5])]
    return return_dict

def up_crunch(dyct,deltaT):
    return_list = []
    for key,value in dyct.iteritems():
        temp,index,count = [],0,0
        if 'A' in value[0][::-1]:
            for sub_item in value[0][::-1]:
                if sub_item == 'A':
                    first = index
                    temp.append(first+1)
                    break
                index+=1
            for sub_item in value[0][::-1][temp[0]-1:]:
                if sub_item == 'A':
                    count +=1
                else:
                    break
            temp.append(count)
            temp.append(deltaT[key])
            return_list.append(temp)
        else:
            temp.append(['NA','NA','NA'])
    return return_list

def down_crunch(dyct,deltaT):
        return_list = []
        for key,value in dyct.iteritems():
            temp,index,count = [],0,0
            nuc = ['C','A','G']
            if any(x in value[1] for x in nuc):
                for sub_item in value[1]:
                    if sub_item in nuc:
                        first = index
                        temp.append(first+1)
                        break
                    index+=1
                for sub_item in value[1][temp[0]-1:]:
                    if sub_item in nuc:
                        count +=1
                    else:
                        break
                temp.append(count)
                temp.append(deltaT[key])
                return_list.append(temp)
            else:
                temp.append(['NA','NA','NA'])
        return return_list

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('file',type=str,help='')
    parser.add_argument('deltaT',type=str,help='')
    args = parser.parse_args()

    seqs,deltaT = seq_strip(args.file),deltaT_dict_build(args.deltaT)
    up_crunched,down_crunched = up_crunch(seqs,deltaT),down_crunch(seqs,deltaT)

    with open('crunched_'+args.file, "wb") as f:
        f.write('dist_up_A,len_up_A,deltaT\n')
        writer = csv.writer(f)
        writer.writerows(up_crunched)

if __name__ == '__main__':
    main()
