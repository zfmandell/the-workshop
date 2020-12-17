#!/usr/bin/env python
from __future__ import division
import sys
from scipy.stats.stats import pearsonr
from itertools import islice
import argparse
from collections import defaultdict
from Bio import SeqIO

def get_covered_transcripts(coverage_fyle):
    '''Reads in a standard overlapped coverage file'''
    info = {}
    with open(coverage_fyle,'r') as f:
        for line in f:
            info[line.strip()] = None
    return info

def filter_dictonary(in_dict,filter_dict):
    '''Removes entries from a dictionary based on coverage'''
    for k in in_dict.keys():
        if k not in filter_dict:
            del in_dict[k]

def threshold_dictonary(in_dict,threshold):
    '''Removes entries from a dictionary based on threshold'''
    """How many zeros you will tolerate"""
    output_dict= {}
    for key,value in in_dict.iteritems():
        if value.count(0.0) <= threshold:
            output_dict[key] = value
    return output_dict

def trend_crunch(list1,list2):
    q = 0

    number_trend = 0

    while q < len(list1):
        one = list1[q]
        two = list2[q]

        if one > 0 and two > 0:
            number_trend += 1
        elif one < 0 and two < 0:
            number_trend += 1

        q += 1

    return number_trend/len(list1)

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

def mod_react(reactvity_dict):
    'outputs -20<x<+10 window for all possible transcripts'

    output_dict = {}

    for key, value in reactvity_dict.iteritems():
        contents = key.split("|")

        new_contents = []

        for item in contents:
            if ":" in str(item):
                new_contents.append(item.split(":"))
            else:
                new_contents.append(item)


        if str(new_contents[1][0]) == 'FPUTR' and int(new_contents[1][1]) == 20:

            output_dict[key] = value[0:31]

        elif str(new_contents[1][0]) == 'FPUTR' and int(new_contents[1][1]) > 20:
            distance = int(new_contents[1][1])-20

            output_dict[key] = value[distance:distance+31]

    return output_dict

def combine_react(ctrl_dict,exp_dict):

    output_dict = {}

    for key_ctrl, value_ctrl in ctrl_dict.iteritems():
        deltas = []
        i = 0
        while i < len(value_ctrl)-1:
            delta = float(exp_dict[key_ctrl][i]) - float(value_ctrl[i])
            deltas.append(delta)
            i += 1
        output_dict[key_ctrl] = deltas

    return output_dict

def binner(threshold_dict):

    scanned = []
    results = []
    tabbed_results =  dict()
    final = []
    final_return = []
    final_return_temp = []
    true_return = []

    for key_one,value_one in threshold_dict.iteritems():
        for key_two,value_two in threshold_dict.iteritems():
            if key_one != key_two:
                to_add = str(key_one)+"+"+str(key_two)

                if to_add not in scanned:
                    scanned.append(to_add)
                    trend_results = (trend_crunch(value_one,value_two))

                    if trend_results > 0.5:
                        results.append(to_add)

    for item in results:
        contents = item.split("+")

        if contents[0] in tabbed_results and contents[1] not in tabbed_results:
            if contents[1] not in tabbed_results[contents[0]]:
                tabbed_results[contents[0]].append(contents[1])
        elif contents[1] in tabbed_results and contents[0] not in tabbed_results:
            if contents[0] not in tabbed_results[contents[1]]:
                tabbed_results[contents[1]].append(contents[0])
        else:
            tabbed_results[contents[0]] = [contents[1]]


    for key,value in tabbed_results.iteritems():
        temp = []
        temp.append(key)
        if value > 1:
            for sub_item in value:
                temp.append(sub_item)
        else:
            temp.append(value[0])

        final.append(temp)

    for item_one in final:
        final[final.index(item_one)] = sorted(item_one)

    for item in final:
        for sub_item in item:
            for second_round in final:
                if sub_item in second_round and final.index(second_round) != final.index(item):
                    for part in item:
                        if part not in final[final.index(second_round)]:
                            final[final.index(second_round)].append(part)
    for item in final:
        final[final.index(item)] = sorted(item)

    for item in final:
        if item not in final_return:
            final_return.append(item)

    for item in final_return:
        temp  = [x.strip(' ') for x in item]
        final_return_temp.append(temp)

    for item in final_return_temp:
        temp = [x.strip("'") for x in item]
        true_return.append(temp)

    return true_return


def main():
    parser = argparse.ArgumentParser(description='Generates changes along n bp from either end of one or more transcripts')
    parser.add_argument('control', type=str, help='Control <.react> file')
    parser.add_argument('experimental',type=str, help='Experimental <.react> file')
    parser.add_argument('-restrict',default = None, help = '<.txt > Limit analysis to these specific transcripts')
    parser.add_argument('-0_threshold',type=int, default=10, help = '<.txt > threshold the number of 0s')
    parser.add_argument('-name',type=str, default='test.txt', help = '<.txt > desired name')

    args = parser.parse_args()

    control,experimental = read_in_derived_reactivities(args.control),read_in_derived_reactivities(args.experimental)

    if args.restrict != None:
        covered = get_covered_transcripts(args.restrict)
        filter_dictonary(control,covered)
        filter_dictonary(experimental,covered)


    SD_control,SD_experimental = mod_react(control),mod_react(experimental)

    SD_delta = combine_react(SD_control,SD_experimental)

    SD_delta_threshold = threshold_dictonary(SD_delta,10)

    binned = binner(SD_delta_threshold)


    with open(args.name,'w') as outp:
        for item in binned:
            outp.write(str(item))
            outp.write("\n")





if __name__ == '__main__':
    main()
