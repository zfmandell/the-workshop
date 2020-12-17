#!/usr/bin/env python

from itertools import islice
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import math
from collections import defaultdict
import argparse


def read_in_derived_reactivities(reactvity_fyle):
    '''Reads in a reactivity file'''
    information = {}
    with open(reactvity_fyle,'r') as f:
        while True:
            next_n_lines = list(islice(f, 2))
            if not next_n_lines:
                break
            transcript,reactivities = [n.strip() for n in next_n_lines]
            information[transcript] = [float(x) if x!='NA' else 0.0 for x in reactivities.split('\t') ]
    return information

def filter_dictonary(in_dict,filter_dict):
    '''Removes entries from a dictionary based on coverage'''
    for k in in_dict.keys():
        if k not in filter_dict:
            del in_dict[k]

def get_covered_transcripts(coverage_fyle):
    '''Reads in a standard overlapped coverage file'''
    info = {}
    with open(coverage_fyle,'r') as f:
        for line in f:
            info[line.strip()] = None
    return info

def subtract_reacts(reacts_A,reacts_B):
    '''subtracts B from A, i.e. A-B'''
    result = []
    spinner = zip(reacts_A,reacts_B)
    for bit in spinner:
        try:
            diff = bit[0]-bit[1]
            result.append(diff)
        except TypeError:
            result.append('NA')
    return result

def get_sum(alist):
    '''returns the sum of a list, removing NAs'''
    return sum([x for x in alist if x != 'NA'])

def get_average(alist):
    return sum(alist)/len(alist)

def scanner(react_list):

    to_return = []
    length = len(react_list)

    q = 0
    while q < length-5:
        summed = get_sum(react_list[q:q+5])
        average = get_average([q,q+4])
        to_add = tuple([summed,average])
        to_return.append(to_add)
        q+=1

    return to_return

def grouper(tuple_list):

    to_return = defaultdict(list)

    for item in tuple_list:
        for sub_item in tuple_list:
            if item[1] != sub_item[1]:
                if int(item[0]) != 0 and int(sub_item[0]) != 0:
                    if abs(int(item[0]) - int(sub_item[0])) < 0.02:
                        to_add = tuple([max([item[1],sub_item[1]]),min([item[1],sub_item[1]])])

                        if max([item[0],sub_item[0]]) not in to_return.keys():
                            to_return[max([item[0],sub_item[0]])] = list(to_add)
                        else:
                            if to_add not in to_return.values():
                                to_return[max([item[0],sub_item[0]])].append(to_add)
    return to_return

def condenser(react_dict):

    to_return = defaultdict(list)

    for key,value in react_dict.iteritems():
        for item in value:
            if isinstance(item, int) != True:
                for sub_item in value:
                    if isinstance(sub_item, int) != True:
                        if item != sub_item:
                            if item[0] in sub_item:
                                print item[0]
                                print sub_item
                                print item
                                to_add = list(item)+list(sub_item)
                                to_add = list(set(to_add)).sort()

                                if key not in to_return.keys():
                                    to_return[key] = list(to_add)
                                else:
                                    to_return[key].append(to_add)


                            elif item[1] in sub_item:
                                to_add = list(item)+list(sub_item)
                                to_add = list(set(to_add)).sort()

                                if key not in to_return.keys():
                                    to_return[key] = list(to_add)
                                else:
                                    to_return[key].append(to_add)

    return to_return


def main():
    parser = argparse.ArgumentParser(description='Generates changes along n bp from either end of one or more transcripts')
    parser.add_argument('control', type=str, help='Control <.react> file')
    parser.add_argument('experimental',type=str, help='Experimental <.react> file')
    parser.add_argument('-restrict',default = None, help = '<.txt > Limit analysis to these specific transcripts')

    args = parser.parse_args()

    control,experimental = read_in_derived_reactivities(args.control),read_in_derived_reactivities(args.experimental)

    if args.restrict != None:
        covered = get_covered_transcripts(args.restrict)
        filter_dictonary(control,covered)
        filter_dictonary(experimental,covered)

    final_dict = {}

    for key,value in control.iteritems():
        difference = subtract_reacts(experimental[key],value)
        scanned = scanner(difference)
        grouped = grouper(scanned)
        condensed = condenser(grouped)
        final_dict[key] = condensed

    print final_dict



if __name__ == '__main__':
    main()
