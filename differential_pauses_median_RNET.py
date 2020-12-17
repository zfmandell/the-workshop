#!/usr/bin/env python3

import argparse
import numpy as np
import collections
import itertools

def genome_yield(fasta_name):
    seq = ''
    with open(fasta_name) as inp:
        for line in inp:
            if line[0] != '>':
                seq = seq+str(line.strip())
    return seq

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases

def pause_place_TU_plus(coord,TU):
    sites = np.array(list(TU.keys()))
    TSS = sites[sites < coord].max()
    return TSS,TU[TSS]

def pause_place_TU_minus(coord,TU):
    sites = np.array(list(TU.keys()))
    TSS = sites[sites > coord].min()
    return TU[TSS],TSS

def pause_place_GFF(gene,strand,GFF3):
    if strand == '+':
        return GFF3[gene][0],GFF3[gene][1]
    else:
        return GFF3[gene][1],GFF3[gene][0]

def read_coverage(cov_fyle):
    #create list of coverage values, position dependent
    return_list = []
    with open(cov_fyle,"r") as inp:
        for line in inp:
            return_list.append(int(line.strip().split('\t')[2]))
    return return_list

def plus_place(pos,lyst):
    for item in sorted(lyst):
        if pos < item:
            return item

def minus_place(pos,lyst):
    for item in sorted(lyst,reverse=True):
        if pos > item:
            return item

def GFF_build(fyle):
    return_dict = {}
    with open(fyle,'r') as inp:
        for line in inp:
            if len(line.strip().split('\t')) > 1:
                if line.strip().split('\t')[2] == 'gene':
                    gene = line.strip().split('\t')[8].split('=')[3].split(';')[0]
                    start = int(line.strip().split('\t')[3])
                    end = int(line.strip().split('\t')[4])
                    strand = line.strip().split('\t')[6]
                    return_dict[gene] = [start,end]
    return return_dict

def TU_build(TSS,TTS):
    TU_plus,TU_minus = {},{}
    TSS_plus,TTS_plus,TSS_minus,TTS_minus = [],[],[],[]
    with open(TTS,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            if line.strip().split(',')[1] == '+':
                TTS_plus.append(int(line.strip().split(',')[0]))
            else:
                TTS_minus.append(int(line.strip().split(',')[0]))
    with open(TSS,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            if line.strip().split(',')[1] == '+':
                TSS_plus.append(int(line.strip().split(',')[0]))
            else:
                TSS_minus.append(int(line.strip().split(',')[0]))

    for item in sorted(TSS_plus):
        TU_plus[item] = plus_place(item,TTS_plus)
    for item in sorted(TSS_minus,reverse=True):
        TU_minus[item] = minus_place(item,TTS_minus)

    empty_plus =  [x[0] for x in [(k, v) for k, v in TU_plus.items() if not v]]
    empty_minus =  [x[0] for x in [(k, v) for k, v in TU_minus.items() if not v]]

    for item in empty_plus:
        TU_plus[item] = sorted(TTS_plus)[0]
    for item in empty_minus:
        TU_minus[item] = sorted(TTS_minus,reverse=True)[0]

    return collections.OrderedDict(sorted(TU_plus.items())),collections.OrderedDict(sorted(TU_minus.items(),reverse=True))

def TPM_scaling(TU_fwd,TU_rev,fwd,rev):
    RPK = []
    for key,value in TU_fwd.items():
        RPK.append(sum(fwd[key:value+1])/(abs(key-value)/1000))
    for key,value in TU_rev.items():
        RPK.append(sum(fwd[value:key+1])/(abs(key-value)/1000))

    return sum(RPK)/1000000

def TPM_calc(coord,strand,TU_fwd,TU_rev,fwd,rev,scaling):

        if strand == '+':
            small,big = pause_place_TU_plus(coord,TU_fwd)
            RPK = sum(fwd[small:big+1])/(abs(small-big)/1000)
            return RPK/scaling
        else:
            small,big = pause_place_TU_minus(coord,TU_rev)
            RPK = sum(rev[small:big+1])/(abs(small-big)/1000)
            return RPK/scaling

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
            return_dict[int(line.strip().split('\t')[0])] = [int(line.strip().split('\t')[1]),float(line.strip().split('\t')[3])]
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
                gene = line.strip().split('\t')[0]
                dist = line.strip().split('\t')[2]
                sense = line.strip().split('\t')[6]
            else:
                peak = int(line.strip().split('\t')[1])
                strand = '-'
                gene = line.strip().split('\t')[0]
                dist = line.strip().split('\t')[2]
                sense = line.strip().split('\t')[6]
            return_dict[peak] = [strand]
            gene_dict[peak] = [gene,sense]

    return return_dict,gene_dict

def pause_rank(pauses_WT,pauses_Mutant,wig_WT,wig_Mutant):
    ranked = {}
    for key,value in pauses_WT.items():
        strand = value[0]
        count_WT = wig_WT[key][0]
        count_Mutant = wig_Mutant[key][0]
        Median_WT = wig_WT[key][1]
        Median_Mutant = wig_Mutant[key][1]
        if Median_WT == 0:
            if strand == '+':
                Median_WT = 1
            else:
                Median_WT = -1
        if abs(Median_Mutant) > 1 and abs(count_Mutant) !=0:
            if np.sign(count_WT) == np.sign(count_Mutant) and np.sign(Median_WT) == np.sign(Median_Mutant):
                Score_WT = round(abs(count_WT)/abs(Median_WT),2)
                Score_Mutant = round(abs(count_Mutant)/abs(Median_Mutant),2)
                ranked[str(key)+'_'+str(strand)] = [abs(count_WT),abs(count_Mutant),abs(Median_WT),abs(Median_Mutant),Score_WT,Score_Mutant]
        elif abs(Median_Mutant) > 1 and abs(count_Mutant) == 0:
            if np.sign(Median_WT) == np.sign(Median_Mutant):
                Score_WT = round(abs(count_WT)/abs(Median_WT),2)
                ranked[str(key)+'_'+str(strand)] = [abs(count_WT),0,abs(Median_WT),abs(Median_Mutant),abs(Score_WT),0]

    for key,value in pauses_Mutant.items():
        strand = value[0]
        count_WT = abs(wig_WT[key][0])
        count_Mutant = abs(wig_Mutant[key][0])
        Median_WT = abs(wig_WT[key][1])
        Median_Mutant = abs(wig_Mutant[key][1])
        if Median_Mutant == 0:
            if strand == '+':
                Median_Mutant = 1
            else:
                Median_Mutant = -1
        if abs(Median_WT) > 1 and abs(count_WT) !=0:
            if np.sign(count_WT) == np.sign(count_Mutant) and np.sign(Median_WT) == np.sign(Median_Mutant):
                Score_WT = round(abs(count_WT)/abs(Median_WT),2)
                Score_Mutant = round(abs(count_Mutant)/abs(Median_Mutant),2)
                ranked[str(key)+'_'+str(strand)] = [abs(count_WT),abs(count_Mutant),abs(Median_WT),abs(Median_Mutant),Score_WT,Score_Mutant]
        elif abs(Median_WT) > 1 and abs(count_WT) == 0:
            if np.sign(Median_WT) == np.sign(Median_Mutant):
                Score_Mutant = round(abs(count_Mutant)/abs(Median_Mutant),2)
                ranked[str(key)+'_'+str(strand)] = [0,abs(count_Mutant),abs(Median_WT),abs(Median_Mutant),0,abs(Score_Mutant)]

    return ranked

def printer(name,dyct,genes_WT,genes_Mutant,WT_fwd,WT_rev,Mutant_fwd,Mutant_rev,GFF3,TU_fwd,TU_rev,scaling_wt,scaling_mut,genome,upstream,downstream):
    with open(name,'w') as outp:
        outp.write('coord,strand,gene,sense,TSS,TTS,start,stop,dist_TSS,dist_TTS,dist_start,dist_stop,WT Count,Mutant Count,WT Median,Mutant Median,WT Score,Mutant Score,WT TPM, Mutant TPM,log2FC,log2(FC)-TPM,log2Diff,sequence,threshold\n')
        for key,value in dyct.items():
            coord = int(key.split('_')[0])
            strand = str(key.split('_')[1])

            if int(coord) in genes_WT.keys():
                gene = str(genes_WT[int(coord)][0])
                sense = str(genes_WT[int(coord)][1])
            else:
                gene = str(genes_Mutant[int(coord)][0])
                sense = str(genes_Mutant[int(coord)][1])

            start,stop = pause_place_GFF(gene,strand,GFF3)

            if strand == '+':
                TSS,TTS = pause_place_TU_plus(coord,TU_fwd)
                dTSS = coord-TSS
                dTTS = coord-TTS
                dstart = coord-start
                dstop = coord-stop
            else:
                TTS,TSS = pause_place_TU_minus(coord,TU_rev)
                dTSS = TSS-coord
                dTTS = TTS-coord
                dstart = start-coord
                dstop = stop-coord

            WT_counts = value[0]
            Mutant_counts = value[1]
            WT_Median = value[2]
            Mutant_Median = value[3]

            WT_Score = value[4]
            if 0 < WT_Score < 1:
                WT_Score == 1

            Mutant_Score = value[5]
            if 0 < Mutant_Score < 1:
                Mutant_score == 1

            WT_TPM = TPM_calc(coord,strand,TU_fwd,TU_rev,WT_fwd,WT_rev,scaling_wt)
            Mutant_TPM = TPM_calc(coord,strand,TU_fwd,TU_rev,Mutant_fwd,Mutant_rev,scaling_mut)

            if WT_TPM > 0:
                WT_Score_TPM = WT_Score/WT_TPM
            else:
                WT_Score_TPM = 0
            if Mutant_TPM > 0:
                Mutant_Score_TPM = Mutant_Score/Mutant_TPM
            else:
                Mutant_Score_TPM = 0

            if WT_Score > 0 and Mutant_Score > 0:
                Log2FC = round(np.log2(Mutant_Score/WT_Score),2)
            else:
                Log2FC = 'NA'
            if WT_Score_TPM > 0 and Mutant_Score_TPM > 0:
                Log2FC_TPM = round(np.log2(Mutant_Score_TPM/WT_Score_TPM),2)
            else:
                Log2FC_TPM = 'NA'

            if abs(WT_Score-Mutant_Score) != 0:
                Log2Diff = round(np.log2(abs(WT_Score-Mutant_Score)),2)
            else:
                Log2Diff = 0
            if abs(WT_Score_TPM-Mutant_Score_TPM) != 0:
                Log2Diff_TPM = round(np.log2(abs(WT_Score_TPM-Mutant_Score_TPM)),2)
            else:
                Log2Diff_TPM = 0

            if Log2FC_TPM != 'NA':
                if abs(Log2FC_TPM) >= 3:
                    differential = 'TRUE'
                else:
                    differential = 'FALSE'

            if int(coord) in genes_WT.keys():
                gene = str(genes_WT[int(coord)][0])
                sense = str(genes_WT[int(coord)][1])
            else:
                gene = str(genes_Mutant[int(coord)][0])
                sense = str(genes_Mutant[int(coord)][1])

            if strand == '+':
                seq = genome[coord-1].upper()
                seq = genome[coord-1-upstream:coord-1].lower()+seq
                seq = seq+genome[coord:coord+downstream].lower()
            else:
                seq = genome[coord-1].upper()
                seq = genome[coord-1-downstream:coord-1].lower()+seq
                seq = seq+genome[coord:coord+upstream].lower()
                seq = reverse_complement(seq)

            if WT_TPM != 0 and Mutant_TPM != 0:
                outp.write(str(coord)+','+strand+','+gene+','+sense+','+str(TSS)+','+str(TTS)+','+str(start)+','+str(stop)+','+str(dTSS)+','+str(dTTS)+','+str(dstart)+','+str(dstop)+\
                ','+str(WT_counts)+','+str(Mutant_counts)+','+str(WT_Median)+','+str(Mutant_Median)+','+str(WT_Score)+','+str(Mutant_Score)+','+str(WT_TPM)+','+str(Mutant_TPM)+\
                ','+str(Log2FC)+','+str(Log2FC_TPM)+','+str(Log2Diff)+','+str(seq)+','+differential+'\n')

def main():
    parser = argparse.ArgumentParser(description='determining pauses that are stimulated/suppressed/independent of mutation')
    parser.add_argument('WT_pauses',type=str,help='Finder 2.1 output from WT [.tab]')
    parser.add_argument('Mutant_pauses',type=str,help='Finder 2.1 output from Mutant [.tab]')
    parser.add_argument('WT_ends',type=str,help='Raw counts at each position in WT [.wig]')
    parser.add_argument('Mutant_ends',type=str,help='Raw counts at each position in Mutant [.wig]')
    parser.add_argument('TSS',type=str,help='List of transcription start sites [.csv]')
    parser.add_argument('TTS',type=str,help='List of transcription termination sites [.csv]')
    parser.add_argument('gene_list',type=str,help='NCBI file denoting coordinates of gene start and stop codons [.gff3]')
    parser.add_argument('cov_fwd_wt',type=str,help='WT fwd strand coverage file from RNET-seq data [.cov]')
    parser.add_argument('cov_rev_wt',type=str,help='WT rev strand coverage file from RNET-seq data [.cov]')
    parser.add_argument('cov_fwd_mut',type=str,help='Mutant fwd strand coverage file from RNET-seq data [.cov]')
    parser.add_argument('cov_rev_mut',type=str,help='Mutant rev strand coverage file from RNET-seq data [.cov]')
    parser.add_argument('genome',type=str,help='genomic reference genome used for mapping [.fasta]')
    parser.add_argument('upstream',type=int,help='# of upstream bases')
    parser.add_argument('downstream',type=int,help='# of downstream bases')
    parser.add_argument('outfile',type=str,help='name of output differential pause stability file')
    args = parser.parse_args()

    #build genome dict
    sequence = genome_yield(args.genome)

    #build pause dict and determine TPM at each 100 bp window around pause site
    pauses_WT,genes_WT = read_pause(args.WT_pauses)
    pauses_Mutant,genes_Mutant = read_pause(args.Mutant_pauses)

    #load up gff file and build TU dict
    GFF3 = GFF_build(args.gene_list)
    TU_fwd,TU_rev = TU_build(args.TSS,args.TTS)

    #clean up TU dicts
    for key,value in dict(TU_fwd).items():
        if value < key:
            del TU_fwd[key]
    for key,value in dict(TU_rev).items():
        if value > key:
            del TU_rev[key]

    #clean up pause dicts
    for key,value in pauses_WT.items():
        if value == '+':
            if key < min(TU_fwd.keys()):
                del pauses_WT[key]
        else:
            if key > max(TU_rev.keys()):
                del pauses_WT[key]

    for key,value in pauses_Mutant.items():
        if value == '+':
            if key < min(TU_fwd.keys()):
                del pauses_Mutant[key]
        else:
            if key > max(TU_rev.keys()):
                del pauses_Mutant[key]

    #load up wig files
    WT_wig,Mutant_wig = read_wig(args.WT_ends),read_wig(args.Mutant_ends)

    #rank all pauses
    ranked = pause_rank(pauses_WT,pauses_Mutant,WT_wig,Mutant_wig)

    #clear large wigs
    WT_wig,Mutant_wig = {},{}

    WT_fwd,WT_rev = read_coverage(args.cov_fwd_wt),read_coverage(args.cov_rev_wt)
    Mutant_fwd,Mutant_rev = read_coverage(args.cov_fwd_mut),read_coverage(args.cov_rev_mut)

    scaling_wt,scaling_mut = TPM_scaling(TU_fwd,TU_rev,WT_fwd,WT_rev),TPM_scaling(TU_fwd,TU_rev,Mutant_fwd,Mutant_rev)

    #print everything out
    printer(args.outfile,ranked,genes_WT,genes_Mutant,WT_fwd,WT_rev,Mutant_fwd,Mutant_rev,GFF3,TU_fwd,TU_rev,scaling_wt,scaling_mut,sequence,args.upstream,args.downstream)

if __name__ == '__main__':
    main()
