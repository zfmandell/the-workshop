#!/usr/bin/env python

from __future__ import division
import sys

def opposite(stryng):
    if stryng == '+':
        return '-'
    else:
        return '+'

def plus_place(pos,lyst):
    for item in lyst:
        if pos < item[0] and item[2] == '+':
            return [item[0],item[3]]

def minus_place(pos,lyst):
    for item in lyst:
        if pos > item[1] and item[2] == '-':
            return [item[1],item[3]]

def normal(pause,TSS,TTS,strand):
    length = 4639675
    if strand == '+':
        if TTS > TSS:
            return round((pause-TSS)/(TTS-TSS)*100)
        else:
            return round((length-TSS+pause)/(length-TSS+TTS)*100)
    else:
        if TTS < TSS:
            return round(((TSS-TTS)-(pause-TTS))/(TSS-TTS)*100)
        else:
            return round((TSS-pause)/(TSS+length-TTS)*100)

def GFF_build(fyle):
    return_list,return_dict = [],{}
    with open(fyle,'r') as inp:
        for line in inp:
            if len(line.strip().split('\t')) > 1:
                if line.strip().split('\t')[2] == 'gene':
                    gene = line.strip().split('\t')[8].split('=')[3].split(';')[0]
                    start = int(line.strip().split('\t')[3])
                    end = int(line.strip().split('\t')[4])
                    strand = line.strip().split('\t')[6]
                    return_list.append([start,end,strand,gene])
                    return_dict[gene] = strand
    return return_list,return_dict

def pause_build(fyle,GFF):
    return_dict = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            if line.strip().split(',')[2] == 'SENSE':
                return_dict[int(line.strip().split(',')[1])] = GFF[line.strip().split(',')[0]]
            else:
                return_dict[int(line.strip().split(',')[1])] = opposite(GFF[line.strip().split(',')[0]])
    return return_dict

def UTR_build(TSS,GFF_plus):
    UTR_plus,UTR_minus = {},{}
    TSS_plus,TSS_minus = [],[]
    GFF_minus = GFF_plus[::-1]
    with open(TSS,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            if line.strip().split(',')[1] == '+':
                TSS_plus.append(int(line.strip().split(',')[0]))
            else:
                TSS_minus.append(int(line.strip().split(',')[0]))

    for item in sorted(TSS_plus):
        UTR_plus[item] = plus_place(item,GFF_plus)

    for item in sorted(TSS_minus,reverse=True):
        UTR_minus[item] = minus_place(item,GFF_minus)

    empty_plus =  [x[0] for x in [(k, v) for k, v in UTR_plus.items() if not v]]
    empty_minus =  [x[0] for x in [(k, v) for k, v in UTR_minus.items() if not v]]

    for item in empty_plus:
        UTR_plus[item] = [GFF_plus[0][0],GFF_plus[0][3]]
    for item in empty_minus:
        UTR_minus[item] = [GFF_minus[0][1],GFF_minus[0][3]]

    return UTR_plus,UTR_minus

def pause_place(pause,plus,minus):
    placement = {}
    for key,value in pause.iteritems():
        if value == '+':
            for sub_key,sub_value in plus.iteritems():
                if sub_key < sub_value[0]:
                    if key > sub_key and key < sub_value[0]:
                        placement[str(key)+'_'+str(sub_key)] = [value,sub_value[0],sub_value[1],key-sub_key,normal(key,sub_key,sub_value[0],'+')]

                else:
                    if key > sub_key and key > sub_value[0]:
                        placement[str(key)+'_'+str(sub_key)] = [value,sub_value[0],sub_value[1],key-sub_key,normal(key,sub_key,sub_value[0],'+')]

        else:
            for sub_key,sub_value in minus.iteritems():
                if sub_key > sub_value[0]:
                    if key < sub_key and key > sub_value[0]:
                        placement[str(key)+'_'+str(sub_key)] = [value,sub_value[0],sub_value[1],sub_key-key,normal(key,sub_key,sub_value[0],'-')]

                else:
                    if key < sub_key and key < sub_value[0]:
                        placement[str(key)+'_'+str(sub_key)] = [value,sub_value[0],sub_value[1],sub_key-key,normal(key,sub_key,sub_value[0],'-')]


    return placement

GFF_coord,GFF = GFF_build(sys.argv[1])
pauses = pause_build(sys.argv[2],GFF)
UTR_plus,UTR_minus = UTR_build(sys.argv[3],GFF_coord)
placed_pauses = pause_place(pauses,UTR_plus,UTR_minus)

with open(sys.argv[4],'w') as outp:
    outp.write('pause,strand,gene,TSS,CDS+1,TSS-dist,UTR_perc\n')
    for key,value in placed_pauses.iteritems():
        outp.write(key.split('_')[0]+','+str(value[0])+','+str(value[2])+','+key.split('_')[1]+','+str(value[1])+','+str(value[3])+','+str(value[4])+'\n')
