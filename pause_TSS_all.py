#!/usr/bin/env python

from __future__ import division
import sys

def opposite(stryng):
    if stryng == '+':
        return '-'
    else:
        return '+'

def plus_place(pos,lyst):
    for item in sorted(lyst):
        if pos < item:
            return item

def minus_place(pos,lyst):
    for item in sorted(lyst,reverse=True):
        if pos > item:
            return item

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

def GFF_dict(fyle):
    return_dict = {}
    with open(fyle,'r') as inp:
        for line in inp:
            if len(line.strip().split('\t')) > 1:
                if line.strip().split('\t')[2] == 'gene':
                    gene = line.strip().split('\t')[8].split('=')[3].split(';')[0]
                    strand = line.strip().split('\t')[6]
                    return_dict[gene] = strand
    return return_dict

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

def TU_build(TSS,TTS):
    TU_plus,TU_minus = {},{}
    TSS_plus,TTS_plus,TSS_minus,TTS_minus = [],[],[],[]
    with open(TTS,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            if line.strip().split(',')[1] == '+':
                TTS_plus.append(int(line.strip().split(',')[0]))
            else:
                TTS_minus.append(int(line.strip().split(',')[0]))
    with open(TSS,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            if line.strip().split(',')[1] == '+':
                TSS_plus.append(int(line.strip().split(',')[0]))
            else:
                TSS_minus.append(int(line.strip().split(',')[0]))

    for item in sorted(TSS_plus):
        TU_plus[item] = plus_place(item,TTS_plus)
    for item in sorted(TSS_minus,reverse=True):
        TU_minus[item] = minus_place(item,TTS_minus)

    empty_plus =  [x[0] for x in [(k, v) for k, v in TU_plus.items() if not v]]
    empty_minus =  [x[0] for x in [(k, v) for k, v in TU_minus.items() if not v]]

    for item in empty_plus:
        TU_plus[item] = sorted(TTS_plus)[0]
    for item in empty_minus:
        TU_minus[item] = sorted(TTS_minus,reverse=True)[0]

    return TU_plus,TU_minus

def pause_place(pause,plus,minus):
    placement = {}
    for key,value in pause.iteritems():
        if value == '+':
            for sub_key,sub_value in plus.iteritems():
                if sub_key < sub_value:
                    if key > sub_key and key < sub_value:
                        placement[str(key)+'_'+str(sub_key)] = [value,sub_value,key-sub_key,normal(key,sub_key,sub_value,'+')]

                else:
                    if key > sub_key and key > sub_value:
                        placement[str(key)+'_'+str(sub_key)] = [value,sub_value,key-sub_key,normal(key,sub_key,sub_value,'+')]

        else:
            for sub_key,sub_value in minus.iteritems():
                if sub_key > sub_value:
                    if key < sub_key and key > sub_value:
                        placement[str(key)+'_'+str(sub_key)] = [value,sub_value,sub_key-key,normal(key,sub_key,sub_value,'-')]
                else:
                    if key < sub_key and key < sub_value:
                        placement[str(key)+'_'+str(sub_key)] = [value,sub_value,sub_key-key,normal(key,sub_key,sub_value,'-')]

    return placement

GFF = GFF_dict(sys.argv[1])
GFF_coord = GFF_dict_coord(sys.argv[1])
pauses = pause_build(sys.argv[2],GFF)
TU_plus,TU_minus = TU_build(sys.argv[3],sys.argv[4])
placed_pauses = pause_place_UTR(pauses,TU_plus,TU_minus)

lost,found = [],[int(x.split('_')[0]) for x in placed_pauses.keys()]
for item in pauses.keys():
    if int(item) not in found:
        placed_pauses[str(item)+'_NA'] = [pauses[item],'NA','NA','NA']

with open(sys.argv[5],'w') as outp:
    outp.write('pause,strand,TSS,TTS,TSS-dist,TU_perc\n')
    for key,value in placed_pauses.iteritems():
        outp.write(key.split('_')[0]+','+str(value[0])+','+key.split('_')[1]+','+str(value[1])+','+str(value[2])+','+str(value[3])+'\n')
