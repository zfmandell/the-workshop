#!/usr/bin/env python

from itertools import islice
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
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

def diff_reacts(reacts_A,reacts_B):
    '''Absolute Value of Difference, two lists'''
    result = []
    spinner = zip(reacts_A,reacts_B)
    for bit in spinner:
        try:
            diff = abs(bit[0]-bit[1])
            result.append(diff)
        except TypeError:
            result.append('NA')
    return result

def get_sum(alist):
    '''returns the sum of a list, removing NAs'''
    if alist:
        return sum([x for x in alist if x != 'NA'])
    else:
        return 0

def find_region(region,tran_name,seq_file,window):
    half_window = window/2

    contents = tran_name.split("|")

    new_contents = []

    for item in contents:
        if ":" in str(item):
            new_contents.append(item.split(":"))
        else:
            new_contents.append(item)

    if str(region).lower() == 'startcodon':

        if str(new_contents[1][0]) == 'FPUTR' and new_contents[1][1] <= half_window:
            return seq_file[:window+3]
        elif str(new_contents[1][0]) == 'FPUTR' and new_contents[1][1] > half_window:
            return seq_file[(int(new_contents[1][1])-half_window):int(new_contents[1][1])+(half_window+3)]
        else:
            return seq_file[:window+3]

    elif str(region).upper() == 'FPUTR':
        if str(new_contents[1][0]) == 'FPUTR':
            return seq_file[:int(new_contents[1][1])]

def hot_stepper(fasta_seqs,control_reacts,experimental_reacts,window=30,region='startcodon'):
    '''Makes lists of reactivity windows'''
    disaster_area = []
    for k, v in fasta_seqs.items():
        try:
            control,experimental = control_reacts[k],experimental_reacts[k]
            change = subtract_reacts(experimental,control)
            abs_change = diff_reacts(experimental,control)
            value_window = find_region(region,k,change,window)
            abs_value_window = find_region(region,k,abs_change,window)
            seq_window = find_region(region,k,v,window)
            if region == 'startcodon':
                tidbit = tuple([get_sum(value_window),get_sum(abs_value_window),seq_window,k])
            elif value_window:
                tidbit = tuple([get_sum(value_window)/len(value_window),get_sum(abs_value_window)/len(abs_value_window),seq_window,k])
            if tidbit not in disaster_area:
                disaster_area.append(tidbit)

        except KeyError:
            continue
    return disaster_area


def ultra_dumper(catastrophe,out='bits.csv'):
    '''Writes out a simple <.csv> file'''
    with open(out,'w') as g:
        g.write('net_change,total_change,seq,transcript\n')
        for item in catastrophe:
            g.write(','.join([str(x) for x in item])+'\n')

def write_out_fasta(info,outfyle='out.fasta',LW=80):
    '''Writes out the <.fasta> file, names are just transcript+step'''
    with open(outfyle,'w') as g:
        for item in info:
            name = item[3]
            seq = item[2]
            g.write('>' + name + '\n')
            for i in xrange(0,len(seq),LW):
                g.write(seq[i:i+LW] + '\n')

def write_out_tab_delim(info,outfyle='out.fasta'):
    with open(outfyle,'w') as g:
        for item in info:
            name = item[3]
            g.write(name+ '\n')

def main():
    #print ''
    #print '\033[1;4;94mStructureFold2:\033[0;0;92m react_windows.py\033[0m'
    #print ''
    parser = argparse.ArgumentParser(description='Generate Windows of changing reactivity, will subtract control from experimental')
    parser.add_argument("control",type=str,help="control <.react> file")
    parser.add_argument("experimental",type=str,help="experimental <.react> file")
    parser.add_argument("fasta",type=str,help="<.fasta> to pull sequences from")
    parser.add_argument('-wlen',type=int, default=30, help='<int> [default = 50] Window Length')
    #parser.add_argument('-wstep',type=int, default=20, help='<int> [default = 20] Window Step')
    parser.add_argument('-outname',type=str,default='windows.csv', help='Change the name of the outfile, overrides default')
    parser.add_argument('-restrict',default = None, help = '<.txt > Limit analysis to these specific transcripts')
    parser.add_argument('-filter_loss',action="store_true",default=False,help = 'Restrict output to largest reactivity losses')
    parser.add_argument('-filter_gain',action="store_true",default=False,help = 'Restrict output to largest reactivity gains')
    parser.add_argument('-filter_delta',action="store_true",default=False,help = 'Restrict output to most changed reactivity')
    parser.add_argument('-perc',type=int, default=25, help='<int> [default = 25] Filter to this percent of windows')
    parser.add_argument('-fastaout',action="store_true",default=False,help = 'Write windows to <.fasta> format as well')
    parser.add_argument('-region',type=str,default='startcodon',help = 'region of interest,options are start codon or FPUTR')
    args = parser.parse_args()
    #
    default_name = '_'.join([args.control.replace('.react',''),args.experimental.replace('.react',''),str(args.wlen)+'win',str(args.region),'FINDER'])+'.csv'
    out_name = default_name if args.outname == 'windows.csv' else args.outname

    #Read in both groups of reactivities, fasta file with sequences
    control_reactivty,experimental_reactivty = read_reactivities(args.control),read_reactivities(args.experimental)
    target_seqs = read_fasta(args.fasta)

    #Truncate seqs to those with good coverage in both conditions if user provides a list
    if args.restrict != None:
        covered = get_covered_transcripts(args.restrict)
        filter_dictonary(target_seqs,covered)

    #Let's get some hot_spots
    loud_noises = hot_stepper(target_seqs,control_reactivty,experimental_reactivty,args.wlen,args.region)

    #Apply filter - they could give multiple, it will work, but it won't really make all that much sense
    if args.filter_loss == True:
        loud_noises = sorted(loud_noises,reverse=False)[0:int((float(args.perc)/100)*len(loud_noises))]
        out_name = out_name.replace('.csv','_'+str(args.perc)+'loss.csv')

    if args.filter_gain == True:
        loud_noises = sorted(loud_noises,reverse=True)[0:int((float(args.perc)/100)*len(loud_noises))]
        out_name = out_name.replace('.csv','_'+str(args.perc)+'gain.csv')

    if args.filter_delta == True:
        loud_noises = sorted(loud_noises,reverse=True,key=lambda x:x[1])[0:int((float(args.perc)/100)*len(loud_noises))]
        out_name = out_name.replace('.csv','_'+str(args.perc)+'delta.csv')

    #Output
    #ultra_dumper(loud_noises,out_name)
    write_out_tab_delim(loud_noises,out_name.replace('.csv','.txt'))
    if args.fastaout == True:
        write_out_fasta(loud_noises,out_name.replace('.csv','.fasta'))

if __name__ == '__main__':
    main()
