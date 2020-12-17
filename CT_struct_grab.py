#!/usr/bin/env python

import glob
import sys
import argparse

def collect_infomation(len_fold):
    #Snags a directory worth of CT files,runs terminal report on each CT file, outputs list of terminal reports
    information = []
    for fyle in glob.glob('*.ct'):
        information.append(terminal_report(fyle,len_fold))
    return information

def terminal_report(fyle,len_fold):
    #generates comprehensive information list of lists, indexed as mentioned above (#1-5)
    return_dict = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        try:
            end,delta = firstline.strip().split()[4].split("_")[0].split(":")[1],firstline.strip().split()[4].split("_")[1].split(":")[1]
            if float(delta) >= 100:
                return_dict[end] = [int(x.split()[4]) for x in [next(inp) for x in xrange(len_fold)]]
        except IndexError:
            end,delta = firstline.strip().split()[1].split("_")[0].split(":")[1],firstline.strip().split()[1].split("_")[1].split(":")[1]
            if float(delta) >= 100:
                return_dict[end] = [0]*len_fold
    return return_dict


def main():
    parser = argparse.ArgumentParser(description='will create .csv file of all delta values at all called peaks ')
    parser.add_argument('len',type=int,help='integer value for length of fold')
    parser.add_argument('outfile',type=str,default=None,help=' base name of output files')
    args = parser.parse_args()

    information = collect_infomation(args.len)
    with open(args.outfile,'w') as outp:
        outp.write("end,")
        for x in xrange(1,args.len+1):
            if x < args.len:
                outp.write("pos_"+str(x)+",")
            elif x == args.len:
                outp.write("pos_"+str(x)+"\n")
        for item in information:
            for key,value in item.iteritems():
                outp.write(str(key)+',')
                for sub_item in value[0:len(value)-1]:
                    outp.write(str(sub_item)+",")
                outp.write(str(value[-1])+"\n")

if __name__ == '__main__':
    main()
