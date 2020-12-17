#!/usr/bin/env python

from itertools import islice
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import math
from collections import defaultdict
import glob
import argparse

def read_fasta(genome_fasta):
    fasta_sequences,fasta_dict =SeqIO.parse(open(genome_fasta),'fasta'),{}
    for fasta in fasta_sequences:
        fasta_dict[fasta.id] = str(fasta.seq)
    return fasta_dict

def read_reactivities(afile):
    '''Reads in a reactivity file'''
    information = {}
    with open(afile,'r') as f:
        while True:
            next_n_lines = list(islice(f, 2))
            if not next_n_lines:
                break
            transcript,reactivities = [n.strip() for n in next_n_lines]
            information[transcript] = [float(x) if x!= 'NA' else 'NA' for x in reactivities.split()]
    return information

def get_covered_transcripts(coverage_fyle):
    '''Reads in a standard overlapped coverage file'''
    info = {}
    with open(coverage_fyle,'r') as f:
        for line in f:
            info[line.strip()] = None
    return info

def filter_dictonary(in_dict,filter_dict):
    '''Removes entries from a dictionary'''
    for k in in_dict.keys():
        if k not in filter_dict:
            del in_dict[k]

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
    return sum([abs(x) for x in alist if x != 'NA'])

def hot_stepper(fasta_seqs,control_reacts,experimental_reacts,window=30,step=30):
    '''Makes lists of reactivity windows'''
    disaster_area = []
    for k, v in fasta_seqs.items():
        try:
            control,experimental = control_reacts[k],experimental_reacts[k]
            change = subtract_reacts(experimental,control)
            value_windows = [change[i:i+window] for i in xrange(0, len(change)-(window-1), step)]
            seq_windows = [v[i:i+window] for i in xrange(0, len(v)-(window-1), step)]
            if len(value_windows) > 1:
                for j in range(0, len(value_windows)):
                    tidbit = tuple([get_sum(value_windows[j]),str(value_windows[j]).replace(",",":"),seq_windows[j],str(j),k])
                    disaster_area.append(tidbit)
        except KeyError:
            continue

    return disaster_area

def binner(noise,ctrl,exp,ctrl_dict,exp_dict,window=30,step=30):
    sums = []
    for item in noise:
        sums.append(item[0])

    top_sum = max(sums)
    bin_length = math.ceil(top_sum)/5

    bins = [x*bin_length for x in range(11)]
    bins = ["%.1f"%item for item in bins]

    for item in noise:
        csvs = glob.glob('./*.csv')
        reacts = glob.glob("./*.react")
        fastas = glob.glob("./*.fasta")
        if isinstance(item[0],float) == True:
            summed = item[0]

        q = 0
        while q < len(bins)-1:
            if summed > float(bins[q]) and summed < float(bins[q+1]):
                if "./"+str(bins[q+1])+"_deltaR_binned.csv" not in csvs:
                    with open(str(bins[q+1])+"_deltaR_binned.csv","w") as outp:
                        outp.write("sum,delta_reactivity,seq,step,transcript\n")
                        for sub_item in item:
                            outp.write(str(sub_item))
                            outp.write(",")
                        outp.write("\n")

                else:
                    with open(str(bins[q+1])+"_deltaR_binned.csv","a") as outp:
                        for sub_item in item:
                            outp.write(str(sub_item))
                            outp.write(",")
                        outp.write("\n")
                if "./"+str(bins[q+1])+"_deltaR_binned.react" not in reacts:
                    with open(str(bins[q+1])+"_deltaR_binned.react","w") as outp:
                        outp.write(str(item[-1]))
                        outp.write("|step:"+str(item[3]))
                        outp.write("\n")
                        to_write = []
                        temp = item[1].replace("[","")
                        temp = temp.replace("]","")
                        temp = temp.replace("'","")
                        temp = temp.replace(" ","")
                        reactivities = temp.split(":")
                        for reactivity in reactivities:
                            if reactivity == 'NA':
                                to_write.append(str(reactivity))
                            else:
                                to_write.append("%.3f"%float(reactivity))

                        outp.write("\t".join(to_write))
                        outp.write("\n")

                else:
                    with open(str(bins[q+1])+"_deltaR_binned.react","a") as outp:
                        outp.write(str(item[-1]))
                        outp.write("|step:"+str(item[3]))
                        outp.write("\n")
                        to_write = []
                        temp = item[1].replace("[","")
                        temp = temp.replace("]","")
                        temp = temp.replace("'","")
                        temp = temp.replace(" ","")
                        reactivities = temp.split(":")
                        for reactivity in reactivities:
                            if reactivity == 'NA':
                                to_write.append(str(reactivity))
                            else:
                                to_write.append("%.3f"%float(reactivity))

                        outp.write("\t".join(to_write))
                        outp.write("\n")


                if "./"+str(bins[q+1])+"_"+str(exp)+"_binned.react" not in reacts:
                    with open(str(bins[q+1])+"_"+str(exp)+"_binned.react","w") as outp:
                        transcript = str(item[-1])
                        step_num = int(item[-2])
                        beginning = step_num*window
                        end = beginning+step

                        out_line_raw = exp_dict[transcript]
                        out_line_filtered = out_line_raw[beginning:end]

                        outp.write(str(item[-1]))
                        outp.write("|step:"+str(item[3]))
                        outp.write("\n")
                        outp.write('\t'.join(map(str,out_line_filtered)))
                        outp.write("\n")


                else:
                    with open(str(bins[q+1])+"_"+str(exp)+"_binned.react","a") as outp:
                        transcript = str(item[-1])
                        step_num = int(item[-2])
                        beginning = step_num*window
                        end = beginning+step

                        out_line_raw = exp_dict[transcript]
                        out_line_filtered = out_line_raw[beginning:end]

                        outp.write(str(item[-1]))
                        outp.write("|step:"+str(item[3]))
                        outp.write("\n")
                        outp.write('\t'.join(map(str,out_line_filtered)))
                        outp.write("\n")

                if "./"+str(bins[q+1])+"_"+str(ctrl)+"_binned.react" not in reacts:
                    with open(str(bins[q+1])+"_"+str(ctrl)+"_binned.react","w") as outp:
                        transcript = str(item[-1])
                        step_num = int(item[-2])
                        beginning = step_num*window
                        end = beginning+step

                        out_line_raw = ctrl_dict[transcript]
                        out_line_filtered = out_line_raw[beginning:end]

                        outp.write(str(item[-1]))
                        outp.write("|step:"+str(item[3]))
                        outp.write("\n")
                        outp.write('\t'.join(map(str,out_line_filtered)))
                        outp.write("\n")


                else:
                    with open(str(bins[q+1])+"_"+str(ctrl)+"_binned.react","a") as outp:
                        transcript = str(item[-1])
                        step_num = int(item[-2])
                        beginning = step_num*window
                        end = beginning+step

                        out_line_raw = ctrl_dict[transcript]
                        out_line_filtered = out_line_raw[beginning:end]

                        outp.write(str(item[-1]))
                        outp.write("|step:"+str(item[3]))
                        outp.write("\n")
                        outp.write('\t'.join(map(str,out_line_filtered)))
                        outp.write("\n")


                if "./"+str(bins[q+1])+"_deltaR_binned.fasta" not in fastas:
                    with open(str(bins[q+1])+"_deltaR_binned.fasta","w") as outp:
                        outp.write(">"+str(item[-1]))
                        outp.write("|step:"+str(item[3]))
                        outp.write("\n")
                        outp.write(str(item[2]))
                        outp.write("\n")


                else:
                    with open(str(bins[q+1])+"_deltaR_binned.fasta","a") as outp:
                        outp.write(">"+str(item[-1]))
                        outp.write("|step:"+str(item[3]))
                        outp.write("\n")
                        outp.write(str(item[2]))
                        outp.write("\n")
            q+= 1


def main():

    parser = argparse.ArgumentParser(description='Generate Windows of changing reactivity, will subtract control from experimental')
    parser.add_argument("control",type=str,help="control <.react> file")
    parser.add_argument("experimental",type=str,help="experimental <.react> file")
    parser.add_argument("fasta",type=str,help="<.fasta> to pull sequences from")
    parser.add_argument("bio_control",type=str,help="biological control condition")
    parser.add_argument("bio_experimental",type=str,help="biological experimental condition")
    parser.add_argument('-wlen',type=int, default=50, help='<int> [default = 30] Window Length')
    parser.add_argument('-restrict',default = None, help = '<.txt > Limit analysis to these specific transcripts')
    parser.add_argument('-perc',type=int, default=50, help='<int> [default = 30] Filter to this percent of windows')
    args = parser.parse_args()

    #Read in both groups of reactivities, fasta file with sequences
    control_reactivty,experimental_reactivty = read_reactivities(args.control),read_reactivities(args.experimental)
    target_seqs = read_fasta(args.fasta)

    ctrl = args.bio_control
    exp = args.bio_experimental

    #Truncate seqs to those with good coverage in both conditions if user provides a list
    if args.restrict != None:
        covered = get_covered_transcripts(args.restrict)
        filter_dictonary(target_seqs,covered)
        filter_dictonary(control_reactivty,covered)
        filter_dictonary(experimental_reactivty,covered)

    #Let's get some hot_spots
    window = int(args.wlen)
    step = int(args.wlen)
    loud_noises = hot_stepper(target_seqs,control_reactivty,experimental_reactivty,window,step)
    loud_noises = sorted(loud_noises,reverse=True)[0:int((float(args.perc)/100)*len(loud_noises))]
    binner(loud_noises,ctrl,exp,control_reactivty,experimental_reactivty,window,step)

if __name__ == '__main__':
    main()
