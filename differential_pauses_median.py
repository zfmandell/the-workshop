#!/usr/bin/env python

from __future__ import division
import argparse
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def median(lst):
    n = len(lst)
    s = sorted(lst)
    return (sum(s[n//2-1:n//2+1])/2.0, s[n//2])[n % 2] if n else None

def read_coverage(cov_fyle):
    #create list of coverage values, position dependent
    return_list = []
    with open(cov_fyle,"r") as inp:
        for line in inp:
            return_list.append(int(line.strip().split('\t')[2]))
    return return_list

def read_wig(wig_fyle):
    #create list of coverage values, position dependent
    return_dict = {}
    with open(wig_fyle,"r") as inp:
        for line in inp:
            return_dict[int(line.strip().split('\t')[0])] = int(line.strip().split('\t')[1])
    return return_dict

def read_pause(fyle):
    return_dict = {}
    gene_dict ={}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            if int(line.strip().split('\t')[8]) > 0:
                peak = int(line.strip().split('\t')[1])
                strand = '+'
                count = abs(int(line.strip().split('\t')[8]))
                gene = line.strip().split('\t')[0]
                dist = line.strip().split('\t')[2]
            else:
                peak = int(line.strip().split('\t')[1])
                strand = '-'
                count = abs(int(line.strip().split('\t')[8]))
                gene = line.strip().split('\t')[0]
                dist = line.strip().split('\t')[2]
            return_dict[peak] = [strand,count]
            gene_dict[peak] = [gene,dist]

    return return_dict,gene_dict

def TPM_build(pauses,fwd,rev):
    return_dict = {}
    for key,value in pauses.iteritems():
        if value[0] == '+':
            med = median(fwd[key-50:key+51])
            return_dict[key] = med
        else:
            med = median(rev[key-50:key+51])
            return_dict[key] = med

    return return_dict

def pause_build(WT_pauses,Mutant_pauses,TPM_WT,TPM_Mutant):
    return_dict = {}
    WT,Mutant = {},{}
    for key,value in WT_pauses.iteritems():
        WT[str(key)+'_'+str(value[0])] = [value[1],value[1]/TPM_WT[key],TPM_WT[key]]
    for key,value in Mutant_pauses.iteritems():
        Mutant[str(key)+'_'+str(value[0])] = [value[1],value[1]/TPM_Mutant[key],TPM_Mutant[key]]
    return_dict['WT'] = WT
    return_dict['Mutant'] = Mutant

    return return_dict

def pause_rank(pauses,floor,fwd_WT,rev_WT,fwd_mutant,rev_mutant,wig_WT,wig_Mutant):
    ranked = {}
    for key,value in pauses['WT'].iteritems():
        if key in pauses['Mutant'].keys():
            ranked[key] = [value[0],pauses['Mutant'][key][0],value[2],pauses['Mutant'][key][2],value[1],pauses['Mutant'][key][1],np.log2(value[1]/pauses['Mutant'][key][1])]
        else:
            pos = int(key.split('_')[0])
            strand = key.split('_')[1]
            WT_Counts = wig_WT[pos]
            Mutant_counts = wig_Mutant[pos]
            WT_TPM = value[2]
            if Mutant_counts != 0:
                if np.sign(WT_Counts) == np.sign(Mutant_counts):
                    if strand == '+':
                        Mutant_TPM = median(fwd_mutant[pos-50:pos+51])
                        if Mutant_TPM > floor:
                            ranked[key] = [abs(WT_Counts),abs(Mutant_counts),abs(value[2]),abs(Mutant_TPM),abs(value[1]),abs(Mutant_counts/Mutant_TPM),np.log2(abs(value[1])/abs(Mutant_counts/Mutant_TPM))]
                    else:
                        Mutant_TPM = median(rev_mutant[pos-50:pos+51])
                        if Mutant_TPM > floor:
                            ranked[key] = [abs(WT_Counts),abs(Mutant_counts),abs(value[2]),abs(Mutant_TPM),abs(value[1]),abs(Mutant_counts/Mutant_TPM),np.log2(abs(value[1])/abs(Mutant_counts/Mutant_TPM))]
                else:
                    if strand == '+':
                        Mutant_TPM = median(fwd_mutant[pos-50:pos+51])
                        if Mutant_TPM > floor:
                            ranked[key] = [abs(WT_Counts),0,abs(value[2]),abs(Mutant_TPM),abs(value[1]),0,'Inf-WT']
                    else:
                        Mutant_TPM = median(rev_mutant[pos-50:pos+51])
                        if Mutant_TPM > floor:
                            ranked[key] = [abs(WT_Counts),0,abs(value[2]),abs(Mutant_TPM),abs(value[1]),0,'Inf-WT']
            else:
                if strand == '+':
                    Mutant_TPM = median(fwd_mutant[pos-50:pos+51])
                    if Mutant_TPM > floor:
                        ranked[key] = [abs(WT_Counts),0,abs(value[2]),abs(Mutant_TPM),abs(value[1]),0,'Inf-WT']
                else:
                    Mutant_TPM = median(rev_mutant[pos-50:pos+51])
                    if Mutant_TPM > floor:
                        ranked[key] = [abs(WT_Counts),0,abs(value[2]),abs(Mutant_TPM),abs(value[1]),0,'Inf-WT']

    for key,value in pauses['Mutant'].iteritems():
        if key not in ranked.keys():
            pos = int(key.split('_')[0])
            strand = key.split('_')[1]
            WT_Counts = wig_WT[pos]
            Mutant_counts = wig_Mutant[pos]
            if WT_Counts != 0:
                if np.sign(WT_Counts) == np.sign(Mutant_counts):
                    if strand == '+':
                        WT_TPM = median(fwd_WT[pos-50:pos+51])
                        if WT_TPM > floor:
                            ranked[key] = [abs(WT_Counts),abs(Mutant_counts),abs(WT_TPM),abs(value[2]),abs(WT_Counts/WT_TPM),abs(value[1]),np.log2(abs(WT_Counts/WT_TPM)/abs(value[1]))]
                    else:
                        WT_TPM = median(rev_WT[pos-50:pos+51])
                        if WT_TPM > floor:
                            ranked[key] = [abs(WT_Counts),abs(Mutant_counts),abs(WT_TPM),abs(value[2]),abs(WT_Counts/WT_TPM),abs(value[1]),np.log2(abs(WT_Counts/WT_TPM)/abs(value[1]))]
                else:
                    if strand == '+':
                        WT_TPM = median(fwd_WT[pos-50:pos+51])
                        if WT_TPM > floor:
                            ranked[key] = [0,abs(Mutant_counts),abs(WT_TPM),abs(value[2]),0,abs(value[1]),'Inf-Mutant']
                    else:
                        WT_TPM = median(rev_WT[pos-50:pos+51])
                        if WT_TPM > floor:
                            ranked[key] = [0,abs(Mutant_counts),abs(WT_TPM),abs(value[2]),0,abs(value[1]),'Inf-Mutant']
            else:
                if strand == '+':
                    WT_TPM = median(fwd_WT[pos-50:pos+51])
                    if WT_TPM > floor:
                        ranked[key] = [0,abs(Mutant_counts),abs(WT_TPM),abs(value[2]),0,abs(value[1]),'Inf-Mutant']
                else:
                    WT_TPM = median(rev_WT[pos-50:pos+51])
                    if WT_TPM > floor:
                        ranked[key] = [0,abs(Mutant_counts),abs(WT_TPM),abs(value[2]),0,abs(value[1]),'Inf-Mutant']

    return ranked

def printer(name,dyct,genes_WT,genes_Mutant):
    with open(name,'w') as outp:
        outp.write('coord,strand,gene,dist,WT_count,Mutant_count,WT_TPM,Mutant_TPM,WT_Score,Mutant_Score,log2FC\n')
        for key,value in dyct.iteritems():
            coord = str(key.split('_')[0])
            strand = str(key.split('_')[1])
            WT_counts = str(value[0])
            Mutant_counts = str(value[1])
            WT_TPM = str(value[2])
            Mutant_TPM = str(value[3])
            WT_Score = str(value[4])
            Mutant_Score = str(value[5])
            Log2FC = str(value[6])
            if int(coord) in genes_WT.keys():
                gene = str(genes_WT[int(coord)][0])
                dist = str(genes_WT[int(coord)][1])
            else:
                gene = str(genes_Mutant[int(coord)][0])
                dist = str(genes_Mutant[int(coord)][1])
            outp.write(coord+','+strand+','+gene+','+dist+','+WT_counts+','+Mutant_counts+','+WT_TPM+','+Mutant_TPM+','+WT_Score+','+Mutant_Score+','+Log2FC+'\n')

def main():
    parser = argparse.ArgumentParser(description='determining pauses that are stimulated/suppressed/independent of mutation')
    parser.add_argument('WT_coverage_fwd',type=str,help='per nucleotide coverage file for WT fwd strand [.cov]')
    parser.add_argument('Mutant_coverage_fwd',type=str,help='per nucleotide coverage file for Mutant fwd strand [.cov]')
    parser.add_argument('WT_coverage_rev',type=str,help='per nucleotide coverage file for WT rev strand [.cov]')
    parser.add_argument('Mutant_coverage_rev',type=str,help='per nucleotide coverage file for Mutant rev strand [.cov]')
    parser.add_argument('WT_pauses',type=str,help='Finder 2.1 output from WT [.tab]')
    parser.add_argument('Mutant_pauses',type=str,help='Finder 2.1 output from Mutant [.tab]')
    parser.add_argument('WT_ends',type=str,help='Raw counts at each position in WT [.wig]')
    parser.add_argument('Mutant_ends',type=str,help='Raw counts at each position in Mutant [.wig]')
    parser.add_argument('protein',type=str,help='protein of interest')
    args = parser.parse_args()

    #load up coverage files
    WT_cov_fwd,Mutant_cov_fwd = read_coverage(args.WT_coverage_fwd),read_coverage(args.Mutant_coverage_fwd)
    WT_cov_rev,Mutant_cov_rev = read_coverage(args.WT_coverage_rev),read_coverage(args.Mutant_coverage_rev)

    #load up wig files
    WT_wig,Mutant_wig = read_wig(args.WT_ends),read_wig(args.Mutant_ends)

    #build pause dict and determine TPM at each 100 bp window around pause site
    pauses_WT,genes_WT = read_pause(args.WT_pauses)
    pauses_Mutant,genes_Mutant = read_pause(args.Mutant_pauses)
    TPM_WT = TPM_build(pauses_WT,WT_cov_fwd,WT_cov_rev)
    TPM_Mutant = TPM_build(pauses_Mutant,Mutant_cov_fwd,Mutant_cov_rev)

    for item in pauses_WT.keys():
        if TPM_WT[item] < 5.0:
            del pauses_WT[item]
            del TPM_WT[item]

    for item in pauses_Mutant.keys():
        if TPM_Mutant[item] < 5.0:
            del pauses_Mutant[item]
            del TPM_Mutant[item]

    #average of lowest TPM value of all pause sites found in WT and Mutant Pause sites
    floor = (min(TPM_WT.values())+min(TPM_Mutant.values()))/2

    print min(TPM_WT.values())
    print min(TPM_Mutant.values())
    print floor

    #build final pause dict
    final_pauses = pause_build(pauses_WT,pauses_Mutant,TPM_WT,TPM_Mutant)

    #rank all pauses
    ranked = pause_rank(final_pauses,floor,WT_cov_fwd,WT_cov_rev,Mutant_cov_fwd,Mutant_cov_rev,WT_wig,Mutant_wig)

    #print everything out
    printer(args.protein+'_ranked_pauses.csv',ranked,genes_WT,genes_Mutant)

if __name__ == '__main__':
    main()
