#!/usr/bin/env python

from glob import glob
import argparse
import collections

def term_build(fyle):
    return_list = []
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            if line.strip().split(",")[1] == 'fwd':
                return_list.append([int(line.strip().split(",")[0]),'+'])
            elif line.strip().split(",")[1] == 'rev':
                return_list.append([int(line.strip().split(",")[0]),'-'])
            else:
                return_list.append([int(line.strip().split(",")[0]),line.strip().split(",")[1]])
    return return_list

def intrinsic_finder(base,query):
    intrinsic = []
    for item in base:
        for sub_item in query:
            if (sub_item[0]-3 <= item[0] <= sub_item[0]+3) and item[1] == sub_item[1]:
                intrinsic.append(sub_item)
    return intrinsic

def perc_term(fyle):
    return_dict= {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            if line.strip().split(",")[1] == 'fwd':
                return_dict[int(line.strip().split(",")[0])] = [float(line.strip().split(",")[2]),'+']
            elif line.strip().split(",")[1] == 'rev':
                return_dict[int(line.strip().split(",")[0])] = [float(line.strip().split(",")[2]),'-']
    return return_dict

def translator(dyct,base):
    return_dict = {}
    for key,value in dyct.iteritems():
        for sub_item in base:
            if sub_item[0]-3 <= key <= sub_item[0]+3 and value[1] == sub_item[1]:
                return_dict[sub_item[0]] = key
    for item in base:
        if item[0] not in return_dict.keys():
            return_dict[item[0]] = 'NA'
    return return_dict

def main():
    parser = argparse.ArgumentParser(description='creates merged index file')
    parser.add_argument('peaks',type=str,default=None,help='name of peak file')
    parser.add_argument('outfile',type=str,help='name of outfile')
    args = parser.parse_args()

    all_ends,all_perc,intrinsic_WT,final_ends,final_perc,translation = {},{},term_build(args.peaks),{},{},{}

    for fyle in glob("*_*.csv"):
        all_ends[fyle.split("_")[0]] = term_build(fyle)
        all_perc[fyle.split("_")[0]] = perc_term(fyle)

    for key,value in all_ends.iteritems():
        final_ends[key] = intrinsic_finder(intrinsic_WT,value)

    for key,value in final_ends.iteritems():
        final_perc[key] = dict([x for x in all_perc[key].items() if x[0] in [y[0] for y in final_ends[key]]])

    WT = collections.OrderedDict(sorted(final_perc['WT'].items()))

    for key,value in final_perc.iteritems():
        print len(value)
        if key != 'WT':
            translation[key] = translator(value,intrinsic_WT)

    with open(args.outfile,'w') as outp:
        outp.write('intrinsic,WT,dA,dG,dR,dAG,dAR,dGR,dAGR\n')
        for item in WT.items():
            outp.write(str(item[0])+','+str(item[1][0])+",")

            if translation['dA'][item[0]] == 'NA':
                outp.write('NA,')
            else:
                outp.write(str(final_perc['dA'][translation['dA'][item[0]]][0])+',')

            if translation['dG'][item[0]] == 'NA':
                outp.write('NA,')
            else:
                outp.write(str(final_perc['dG'][translation['dG'][item[0]]][0])+',')

            if translation['dR'][item[0]] == 'NA':
                outp.write('NA,')
            else:
                outp.write(str(final_perc['dR'][translation['dR'][item[0]]][0])+',')

            if translation['dAG'][item[0]] == 'NA':
                outp.write('NA,')
            else:
                outp.write(str(final_perc['dAG'][translation['dAG'][item[0]]][0])+',')

            if translation['dAR'][item[0]] == 'NA':
                outp.write('NA,')
            else:
                outp.write(str(final_perc['dAR'][translation['dAR'][item[0]]][0])+',')

            if translation['dGR'][item[0]] == 'NA':
                outp.write('NA,')
            else:
                outp.write(str(final_perc['dGR'][translation['dGR'][item[0]]][0])+',')

            if translation['dAGR'][item[0]] == 'NA':
                outp.write('NA\n')
            else:
                outp.write(str(final_perc['dAGR'][translation['dAGR'][item[0]]][0])+'\n')


if __name__ == '__main__':
    main()
