#!/usr/bin/env python

import argparse
from collections import defaultdict


def csv_to_dict(csv):
    to_return = defaultdict(list)
    to_return_sorted = {}
    with open(csv,'rU') as f:
        firstline = f.readline()
        for line in f:
            line = line.strip()
            contents = line.split(",")
            transcript = contents[3]
            step = contents[4]
            value = contents[1]

            if transcript not in to_return.keys():
                to_return[transcript] = list(step)
                to_return[transcript].append(value)
            else:
                to_return[transcript].append(step)
                to_return[transcript].append(value)

    return to_return

def compare_print(dict_list,name):
    with open(name,'w') as outp:
        outp.write("transcript,step,WT,deaD,cspA\n")

        for key,value in dict_list[0].iteritems():
            outp.write(str(key)+",")
            q = 1
            value_WT = value
            value_deaD = dict_list[1][key]
            value_cspA = dict_list[2][key]
            for item in value_WT:
                if q&1:
                    outp.write(str(item)+",")
                    try:
                        index_WT = value_WT.index(item)
                        total_WT = value_WT[index_WT+1]
                    except ValueError:
                        total_WT = 'NA'

                    try:
                        index_deaD = value_deaD.index(item)
                        total_deaD = value_deaD[index_deaD+1]
                    except ValueError:
                        total_deaD = 'NA'

                    try:
                        index_cspA = value_cspA.index(item)
                        total_cspA = value_deaD[index_cspA+1]
                    except ValueError:
                        total_cspA = 'NA'


                    outp.write(total_WT)
                    outp.write(",")
                    outp.write(total_deaD)
                    outp.write(",")
                    outp.write(total_cspA)
                    outp.write("\n")

                    q += 1

def main():
    parser = argparse.ArgumentParser(description='Generate merged overlap file')
    parser.add_argument("WT",type=str,help="overlap file 1")
    parser.add_argument("deaD",type=str,help="overlap file 2")
    parser.add_argument("cspA",type=str,help="overlap file 2")
    parser.add_argument('name',type=str, help='name of new overlap file')

    args = parser.parse_args()

    WT_dict = csv_to_dict(args.WT)
    deaD_dict = csv_to_dict(args.deaD)
    cspA_dict = csv_to_dict(args.cspA)

    thermo_dicts = [WT_dict,deaD_dict,cspA_dict]

    compare_print(thermo_dicts,args.name)


if __name__ == '__main__':
    main()
