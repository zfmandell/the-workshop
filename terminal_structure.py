#!/usr/bin/env python

"""
This Script takes a CT folder from batch_fold_rna.py module of StructureFold2
outputs basic stats about all structure proximal to 3' end of folder
eg)
1) distance of nearest Structure (will discount base pairing in 5' terminal 50 nt, due to dubious nature of this interaction)
2) length of base-pairing for nearest structure (ignores one nucleotide mismatch)
3) how structured is the 3' terminal 50 nt (percent base paired,will discount base pairing in 5' terminal 50 nt, due to dubious nature of this interaction)
and basic stats on all of these, mean, median, mode, std dev
lastly, will rank each terminal end as a function of item #5, will place each end in genomic context
"""

from __future__ import division
import os
import glob
import operator
import numpy as np
from scipy import stats
import argparse

def read_gff3(gff3):
    #reads in fasta, returns dict of sequences, dict of strand info for gene
    return_dict = {}
    with open(gff3,'r') as inp:
        for line in inp:
            if line[0] != '#' and len(line.strip().split()) != 0:
                if str(line.strip().split()[2]) == 'gene':
                    return_dict[(line.strip().split()[8].split(";")[2].split("=")[1]).replace("-",".")] = [int(line.strip().split()[3]),int(line.strip().split()[4])]

    return sorted(return_dict.items(), key=operator.itemgetter(1))

def collect_infomation():
    #Snags a directory worth of CT files,runs terminal report on each CT file, outputs list of terminal reports
    information = []
    for fyle in glob.glob('*.ct'):
        information.append(terminal_report(fyle))
    return information

def gene_compile(ends,seqs):
    #returns dict where coordinate : genomical_feature (either CDS, or between two CDS)
    feature_dict = {}

    for item in ends:
        for seq in seqs:
            if item < seqs[0][1][0]:
                feature_dict[str(item)] = "-".join([seqs[0][0],seqs[-1][0]])
            elif seq[1][0] <= item <= seq[1][1]:
                feature_dict[str(item)] = str(seq[0])
            elif seqs.index(seq) < len(seqs) - 1:
                if seq[1][1] < item < seqs[seqs.index(seq)+1][1][0]:
                    feature_dict[str(item)] = "-".join([seq[0],seqs[seqs.index(seq)+1][0]])
            elif seqs.index(seq) == len(seqs) - 1:
                if seq[1][1] < item:
                    feature_dict[str(item)] = "-".join([seq[0],seqs[0][0]])

    return feature_dict

def structured_context(info_reports):
    #takes information reports from collect information, returns dict where terminal_structure (item 3) : coordinate of 3' ends
    return_dict = {}

    for item in info_reports:
        return_dict[item[0]] = item[3]

    return return_dict

def terminal_report(fyle):
    #generates comprehensive information list of lists, indexed as mentioned above (#1-5)
    struct = []

    with open(fyle,'r') as inp:
        firstline = inp.readline()
        end = firstline.strip().split()[1].split("_")[0].split(":")[1]

        for line in inp:
            struct.append(line.strip().split()[4])

    rev_struct = list(reversed(struct))

    q = 0
    for item in rev_struct:
        q+=1
        if int(item) > 51:
            distance_start = q
            break

    i = 0
    it = iter(rev_struct[distance_start-1:])
    for x in it:
        i+=1
        if int(x) == 0 and int(next(it)) == 0:
            distance_end = q+i
            break

    distance = distance_end-distance_start

    struct_tab = 0
    for item in rev_struct[0:51]:
        if int(item) > 51:
            struct_tab+=1
    struct = struct_tab/50

    return [end,distance_start,distance,struct]

def stats(info_sublist):
    #generates basic stats when passed a list, returns list of [mean,min,max,median,std_dev] of that list
    arrayed = np.asarray(info_sublist)
    return [np.mean(arrayed),np.amin(arrayed),np.amax(arrayed),np.median(arrayed),np.std(arrayed)]

def writer(information,outfile,features,end_structure):

    dist = stats([int(item[1]) for item in information])
    length_stem = stats([item[2] for item in information])
    terminal_structure = stats([item[3] for item in information])

    ordered_ends = sorted(end_structure.items(), key=operator.itemgetter(1),reverse=True)

    with open(outfile,'w') as outp:
        outp.write("distance of nearest Structure (discounting base pairing in 5' terminal 50 nt, due to dubious nature of this interaction)\n")
        outp.write("mean,min,max,median,std_dev\n")
        for item in dist[0:-2]:
            outp.write(str(item)+",")
        outp.write(str(dist[-1]))
        outp.write("\n")
        outp.write("\n")

        outp.write("length of base-pairing for nearest structure (ignores one nucleotide mismatch)\n")
        outp.write("mean,min,max,median,std_dev\n")
        for item in length_stem[0:-2]:
            outp.write(str(item)+",")
        outp.write(str(length_stem[-1]))
        outp.write("\n")
        outp.write("\n")

        outp.write("how structured is the 3' terminal 50 nt (percent base paired,will discount base pairing in 5' terminal 50 nt, due to dubious nature of this interaction)\n")
        outp.write("mean,min,max,median,std_dev\n")
        for item in terminal_structure[0:-2]:
            outp.write(str(item)+",")
        outp.write(str(terminal_structure[-1]))
        outp.write("\n")
        outp.write("\n")

        q = 1
        outp.write("3' terminal 50 nt ordered as a function of structure (top = most structured, bottom = least structured )\n")
        outp.write("rank,percent_stranded,coordinate,gene_feature\n")
        for item in ordered_ends:
            outp.write(str(q)+','+str(item[1])+','+str(item[0])+','+str(features[item[0]])+"\n")
            q+=1

def main():
    parser = argparse.ArgumentParser(description='will create .txt file summarizing CT file created from batch_fold_rna.py')
    parser.add_argument('genbank',type=str,help='<.gff3> reference genome from genbank in gff3 format')
    parser.add_argument('outfile',type=str,help='output file name')
    args = parser.parse_args()

    information = collect_infomation()
    ends = [float(item[0]) for item in information]
    seqs = read_gff3(args.genbank)
    features = gene_compile(ends,seqs)
    end_stucture = structured_context(information)
    writer(information,args.outfile,features,end_stucture)

if __name__ == '__main__':
    main()
