#!/usr/bin/env python

"""
calculates the delta value at each position of a wig file eg) [x bp upstream] - [x bp downstream]. outputs delta wig file

x set by user, default is 15 bp
"""

import argparse
import glob
from os.path import basename

def read_coverage(wig_fyle):
    #create list of coverage values, position dependent
    with open(wig_fyle,"r") as inp:
        firstline = inp.readline()
        output_list = [[float(i) for ind,i in enumerate(line.strip().split("\t")) if ind == 3] for line in inp.readlines() if str(line[0]) != '#']
        flat_list = [item for sublist in output_list for item in sublist]
    print output
    return flat_list

def delta_calc(coverage_list,windows):
    #creates list of delta value, position dependent, first and third conditional is due to circular nature of bacterial chromosome
    #'upstream' of the 5' end of the list is the end of the list (reverse)
    #'upstream' of the 3' end of the list is the beginning of the list
    output = []
    counter = 0
    length = len(coverage_list)

    while counter < length:
        if counter < windows:
            gap = windows - counter
            gapped = sum(coverage_list[-gap:])
            output.append((gapped+sum(coverage_list[:counter]))-sum(coverage_list[counter:counter+windows]))
            counter += 1
        elif counter >= windows and counter <= (length - windows):
            output.append(sum(coverage_list[counter-windows:counter])-sum(coverage_list[counter:counter+windows]))
            counter += 1
        elif counter > (length - windows):
            gap = windows - (length-counter)
            gapped = sum(coverage_list[:gap])
            output.append(sum(coverage_list[counter-windows:counter])-(gapped+sum(coverage_list[counter:])))
            counter += 1
    print output
    return output

def writer(wig_fyle,delta_list):
    #creates delta wig file, based off values in wig file
    new_name = 'delta_'+str(wig_fyle[:-4])+".bedgraph"
    counter = 0
    with open(wig_fyle,'r') as inp:
        firstline = inp.readline()
        with open(new_name,"w") as outp:
            outp.write(firstline)
            while counter < len(delta_list)-1:
                for line in inp:
                    if str(line[0]) == '#':
                        outp.write(line)
                    else:
                        line = line.strip().split("\t")
                        outp.write("\t".join(line[0:3]))
                        outp.write("\t")
                        if delta_list[counter] >= 0:
                            outp.write(str(delta_list[counter]))
                        elif delta_list[counter] < 0:
                            outp.write('0')
                        outp.write("\n")
                        counter += 1


def main():
    parser = argparse.ArgumentParser(description='calculates the delta value at each position [ x upstream - x downstream].')
    parser.add_argument('-window_size',type=int,default=15,help='File 1, file you wish to supplement')
    parser.add_argument('-file', default = None, help='Specific <.wig>s to use', nargs='+',dest='wigs')
    parser.add_argument('-batchdir',action="store_true",default=False,help = 'Use all coverage <.wig> in the directory')
    args = parser.parse_args()

    if args.file != None:
        for fyle in args.file:
            print fyle
            sub_data = read_coverage(fyle)
            sub_delta = delta_calc(sub_data,args.window_size)
            print len(sub_delta)
            print len(sub_data)
            writer(fyle,sub_delta)

    if args.batchdir == True:
        for fyle in sorted(glob.glob('*.wig')):
            sub_data = read_coverage(fyle)
            sub_delta = delta_calc(sub_data,args.window_size)
            writer(fyle,sub_delta)


if __name__ == '__main__':
    main()
