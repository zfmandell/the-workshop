#!/usr/bin/env python

import argparse
import glob
import operator
import sys


def translator(lyst,base):
    for item in base:
        if item[0]-3 <= lyst[0] <= item[0]+3 and lyst[1] == item[1][0]:
            return item[0]

final = {}
for fyle in glob.glob('./*rep*'):
    to_add = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            to_add[int(line.strip().split(",")[0])] = [line.strip().split(",")[1],float(line.strip().split(",")[2]),float(line.strip().split(",")[3]),float(line.strip().split(",")[4])]
    final[fyle] = to_add

translated = {}
for key,value in final.iteritems():
    to_add = {}
    if key != './WT_term_eff_rep.csv':
        for sub_key,sub_value in value.iteritems():
            to_add[sub_key] = translator([sub_key,sub_value[0]],final['./WT_term_eff_rep.csv'].items())
    translated[key] = to_add

final_translated = {}
for key,value in final.iteritems():
    if key != './WT_term_eff_rep.csv':
        to_add = {}
        for sub_key,sub_value in value.iteritems():
            if translated[key][sub_key] != None:
                to_add[translated[key][sub_key]] = sub_value
            else:
                to_add[sub_key] = sub_value
        final_translated[key] = to_add
    else:
        final_translated[key] = value

ends = []
for key,value in final_translated.iteritems():
    if key == './WT_term_eff_rep.csv':
        for sub_key,sub_value in value.iteritems():
            if [sub_key,sub_value[0]] not in ends:
                ends.append([sub_key,sub_value[0]])
    elif key == './dA_term_eff_rep.csv':
        for sub_key,sub_value in value.iteritems():
            if [sub_key,sub_value[0]] not in ends:
                ends.append([sub_key,sub_value[0]])
    elif key == './dG_term_eff_rep.csv':
        for sub_key,sub_value in value.iteritems():
            if [sub_key,sub_value[0]] not in ends:
                ends.append([sub_key,sub_value[0]])
    elif key == './dAG_term_eff_rep.csv':
        for sub_key,sub_value in value.iteritems():
            if [sub_key,sub_value[0]] not in ends:
                ends.append([sub_key,sub_value[0]])
print len(ends)

with open(sys.argv[1],'w') as outp:
    outp.write('POT,WT_1,WT_2,WT_3,dA_1,dA_2,dA_3,dG_1,dG_2,dG_3,dAG_1,dAG_2,dAG_3\n')
    for item in sorted(ends):
        outp.write(str(item[0]))
        outp.write(",")
        if item[0] in final_translated['./WT_term_eff_rep.csv'].keys() and final_translated['./WT_term_eff_rep.csv'][item[0]][0] == item[1]:
            outp.write(str(final_translated['./WT_term_eff_rep.csv'][item[0]][1]))
            outp.write(",")
            outp.write(str(final_translated['./WT_term_eff_rep.csv'][item[0]][2]))
            outp.write(",")
            outp.write(str(final_translated['./WT_term_eff_rep.csv'][item[0]][3]))
            outp.write(",")
        else:
            outp.write('0,0,0,')
        if item[0] in final_translated['./dA_term_eff_rep.csv'].keys() and final_translated['./dA_term_eff_rep.csv'][item[0]][0] == item[1]:
            outp.write(str(final_translated['./dA_term_eff_rep.csv'][item[0]][1]))
            outp.write(",")
            outp.write(str(final_translated['./dA_term_eff_rep.csv'][item[0]][2]))
            outp.write(",")
            outp.write(str(final_translated['./dA_term_eff_rep.csv'][item[0]][3]))
            outp.write(",")
        else:
            outp.write('0,0,0,')
        if item[0] in final_translated['./dG_term_eff_rep.csv'].keys() and final_translated['./dG_term_eff_rep.csv'][item[0]][0] == item[1]:
            outp.write(str(final_translated['./dG_term_eff_rep.csv'][item[0]][1]))
            outp.write(",")
            outp.write(str(final_translated['./dG_term_eff_rep.csv'][item[0]][2]))
            outp.write(",")
            outp.write(str(final_translated['./dG_term_eff_rep.csv'][item[0]][3]))
            outp.write(",")
        else:
            outp.write('0,0,0,')
        if item[0]in final_translated['./dAG_term_eff_rep.csv'].keys() and final_translated['./dAG_term_eff_rep.csv'][item[0]][0] == item[1]:
            outp.write(str(final_translated['./dAG_term_eff_rep.csv'][item[0]][1]))
            outp.write(",")
            outp.write(str(final_translated['./dAG_term_eff_rep.csv'][item[0]][2]))
            outp.write(",")
            outp.write(str(final_translated['./dAG_term_eff_rep.csv'][item[0]][3]))
            outp.write("\n")
        else:
            outp.write('0,0,0\n')

        """if item[0] in final_translated['./dR_term_eff_rep.csv'].keys() and final_translated['./dR_term_eff_rep.csv'][item[0]][0] == item[1]:
            outp.write(str(final_translated['./dR_term_eff_rep.csv'][item[0]][1]))
            outp.write(",")
            outp.write(str(final_translated['./dR_term_eff_rep.csv'][item[0]][2]))
            outp.write(",")
            outp.write(str(final_translated['./dR_term_eff_rep.csv'][item[0]][3]))
            outp.write("\n")
        else:
            outp.write('0,0,0\n')
        outp.write(str(item[0]))
        outp.write(",")
        outp.write(item[1])
        outp.write(",")
        outp.write('dAR,')
        if item[0] in final_translated['./dAR_term_eff_rep.csv'].keys() and final_translated['./dAR_term_eff_rep.csv'][item[0]][0] == item[1]:
            outp.write(str(final_translated['./dAR_term_eff_rep.csv'][item[0]][1]))
            outp.write(",")
            outp.write(str(final_translated['./dAR_term_eff_rep.csv'][item[0]][2]))
            outp.write(",")
            outp.write(str(final_translated['./dAR_term_eff_rep.csv'][item[0]][3]))
            outp.write("\n")
        else:
            outp.write('0,0,0\n')
        outp.write(str(item[0]))
        outp.write(",")
        outp.write(item[1])
        outp.write(",")
        outp.write('dGR,')
        if item[0] in final_translated['./dGR_term_eff_rep.csv'].keys() and final_translated['./dGR_term_eff_rep.csv'][item[0]][0] == item[1]:
            outp.write(str(final_translated['./dGR_term_eff_rep.csv'][item[0]][1]))
            outp.write(",")
            outp.write(str(final_translated['./dGR_term_eff_rep.csv'][item[0]][2]))
            outp.write(",")
            outp.write(str(final_translated['./dGR_term_eff_rep.csv'][item[0]][3]))
            outp.write("\n")
        else:
            outp.write('0,0,0\n')
        outp.write(str(item[0]))
        outp.write(",")
        outp.write(item[1])
        outp.write(",")
        outp.write('dAG,')
        if item[0]in final_translated['./dAG_term_eff_rep.csv'].keys() and final_translated['./dAG_term_eff_rep.csv'][item[0]][0] == item[1]:
            outp.write(str(final_translated['./dAG_term_eff_rep.csv'][item[0]][1]))
            outp.write(",")
            outp.write(str(final_translated['./dAG_term_eff_rep.csv'][item[0]][2]))
            outp.write(",")
            outp.write(str(final_translated['./dAG_term_eff_rep.csv'][item[0]][3]))
            outp.write("\n")
        else:
            outp.write('0,0,0\n')
        outp.write(str(item[0]))
        outp.write(",")
        outp.write(item[1])
        outp.write(",")
        outp.write('dAGR,')
        if item[0] in final_translated['./dAGR_term_eff_rep.csv'].keys() and final_translated['./dAGR_term_eff_rep.csv'][item[0]][0] == item[1]:
            outp.write(str(final_translated['./dAGR_term_eff_rep.csv'][item[0]][1]))
            outp.write(",")
            outp.write(str(final_translated['./dAGR_term_eff_rep.csv'][item[0]][2]))
            outp.write(",")
            outp.write(str(final_translated['./dAGR_term_eff_rep.csv'][item[0]][3]))
            outp.write("\n")
        else:
            outp.write('0,0,0\n')"""
