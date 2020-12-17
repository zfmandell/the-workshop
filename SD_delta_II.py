#!/usr/bin/env python

from __future__ import division
from itertools import islice
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import argparse

def read_fasta(genome_fasta):
    fasta_sequences,fasta_dict =SeqIO.parse(open(genome_fasta),'fasta'),{}
    for fasta in fasta_sequences:
        fasta_dict[str(fasta.id)] = str(fasta.seq)
    return fasta_dict

def filter_dictonary(in_dict,filter_dict):
    '''Removes entries from a dictionary based on coverage'''
    for k in in_dict.keys():
        if k not in filter_dict:
            del in_dict[k]

def get_covered_transcripts(coverage_fyle):
    '''Reads in a non standard overlapped coverage file from SD_delta'''
    temp = []
    info = {}
    with open(coverage_fyle,'r') as f:
        for line in f:
            line = line.strip()
            line = line.split(",")
            temp.append(line)

        strongest = max(temp,key=len)

        return_temp  = [x.strip(' ') for x in strongest]
        return_temp  = [x.strip("'") for x in return_temp]

        for item in return_temp:
            info[item] = None
    return info

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
        exp_value = exp_dict[key_ctrl]
        deltas = []
        i = 0
        while i < len(value_ctrl)-1:
            delta = float(exp_value[i]) - float(value_ctrl[i])
            deltas.append(delta)
            i += 1
        output_dict[key_ctrl] = deltas

    return output_dict

def averaged_delta(ctrl_dict,exp_dict,output_fyle):

    output_list_ctrl = []
    output_list_exp = []
    merged = []

    q = 0

    values_ctrl = ctrl_dict.values()
    values_exp = exp_dict.values()

    length = len(values_ctrl[0])

    while q < length:

        temp_ctrl  = [item_ctrl[q] for item_ctrl in values_ctrl]
        temp_exp  = [item_exp[q] for item_exp in values_exp]


        temp_ctrl_average = sum(temp_ctrl)/float(len(temp_ctrl))
        temp_exp_average = sum(temp_exp)/float(len(temp_exp))

        temp_merged = [temp_ctrl_average,temp_exp_average]

        merged.append(temp_merged)

        q += 1


    i = -20

    with open(output_fyle,"w") as outp:
        outp.write("position,averaged_reactivity_ctrl,averaged_reactivity_exp\n")

        for item in merged:
            outp.write(str(i)+","+str(item[0])+","+str(item[1]))
            outp.write("\n")

            i += 1

def folding_prep_react(react_fyle):

    with open(react_fyle,"r") as inp:

        output_name = 'SD_'+str(react_fyle)
        with open(output_name,"a+") as outp:

            integer = 1

            prevline = ''

            for line in inp:
                line = line.strip()

                if integer % 2 != 0:
                    prevline = line
                    integer += 1

                elif integer % 2 == 0:

                    contents = prevline.split("|")

                    new_contents = []

                    for item in contents:
                        if ":" in str(item):
                            new_contents.append(item.split(":"))
                        else:
                            new_contents.append(item)

                    if str(new_contents[1][0]) == 'FPUTR' and int(new_contents[1][1]) == 20:
                        base_contents = line.split("\t")

                        new_line = base_contents[:31]

                        header =  "|".join(contents)
                        outp.write(header)
                        outp.write('\n')
                        outp.write('\t'.join(new_line))
                        outp.write('\n')

                    elif str(new_contents[1][0]) == 'FPUTR' and int(new_contents[1][1]) > 20:
                        base_contents = line.split("\t")

                        distance = int(new_contents[1][1])-20

                        new_line = base_contents[distance:distance+31]

                        header =  "|".join(contents)
                        outp.write(header)
                        outp.write('\n')
                        outp.write('\t'.join(new_line))
                        outp.write('\n')

                    integer +=1

def folding_prep_fasta(fasta_fyle,name,coverage):
    for record in SeqIO.parse(fasta_fyle, "fasta"):

        if str(record.id) in coverage.keys():

            output_name = str(name)+str(fasta_fyle)

            contents = str(record.id).split("|")

            new_contents = []

            for item in contents:
                if ":" in str(item):
                    new_contents.append(item.split(":"))
                else:
                    new_contents.append(item)

            if str(new_contents[1][0]) == 'FPUTR' and int(new_contents[1][1]) == 20:

                sequence = str(record.seq)[0:31]

                with open(output_name,'a') as fasta:
                    fasta.write(">"+str(record.id))
                    fasta.write("\n")
                    fasta.write(sequence)
                    fasta.write("\n")

            elif str(new_contents[1][0]) == 'FPUTR' and int(new_contents[1][1]) > 20:

                distance = int(new_contents[1][1])-20

                sequence = str(record.seq)[distance:distance+31]

                with open(output_name,'a+') as fasta:
                    fasta.write(">"+str(record.id))
                    fasta.write("\n")
                    fasta.write(sequence)
                    fasta.write("\n")


def main():
    parser = argparse.ArgumentParser(description='processing from SD_delta.py')
    parser.add_argument('control', type=str, help='Control <.react> file')
    parser.add_argument('experimental',type=str, help='Experimental <.react> file')
    parser.add_argument('-restrict',default = None, help = '<.txt > Limit analysis to these specific transcripts')
    parser.add_argument('-name',type=str, default='test.txt', help = '<.txt > desired name')
    parser.add_argument('-fasta',type=str, default=None, help = '<.txt > desired fasta')

    args = parser.parse_args()

    control,experimental = read_in_derived_reactivities(args.control),read_in_derived_reactivities(args.experimental)

    covered = get_covered_transcripts(args.restrict)
    filter_dictonary(control,covered)
    filter_dictonary(experimental,covered)

    SD_control,SD_experimental = mod_react(control),mod_react(experimental)

    SD_delta = combine_react(SD_control,SD_experimental)

    averaged_name = str(args.name)+".csv"

    averaged_delta(SD_control,SD_experimental,averaged_name)

    folding_prep_react(args.control)
    folding_prep_react(args.experimental)
    folding_prep_fasta(args.fasta,args.name,covered)


if __name__ == '__main__':
    main()
