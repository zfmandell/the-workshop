#!/usr/bin/env python3

from __future__ import division
import argparse
import numpy as np
import glob
import operator
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from bioinfokit import analys, visuz
import pandas as pd


def median(lst):
    n = len(lst)
    s = sorted(lst)
    return (sum(s[n//2-1:n//2+1])/2.0, s[n//2])[n % 2] if n else None

def read_wig(wig_fyle):
    #create list of coverage values, position dependent
    return_dict = {}
    with open(wig_fyle,"r") as inp:
        for line in inp:
            return_dict[int(line.strip().split('\t')[0])] = float(line.strip().split('\t')[2])
    return return_dict

def read_tab(tab_fyle):
    return_dict = {}
    gene_dict ={}
    with open(tab_fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            if int(line.strip().split('\t')[8]) > 0:
                peak = int(line.strip().split('\t')[1])
                strand = '+'
                gene = line.strip().split('\t')[0]
                dist = line.strip().split('\t')[2]
            else:
                peak = int(line.strip().split('\t')[1])
                strand = '-'
                gene = line.strip().split('\t')[0]
                dist = line.strip().split('\t')[2]
            return_dict[peak] = strand
            gene_dict[peak] = [gene,dist]

    return return_dict,gene_dict

def pause_rank(pauses,wt_wig,mutant_wig):
    ranked_log2,ranked_pval,raw_val = {},{},{}
    wt_values,mutant_values = {},{}
    for key,value in pauses.items():
        strand_pause = value
        hashed = str(key)+'_'+str(value)
        for item in wt_wig:
            score = item[key]

            if score == 0 and strand_pause == "+":
                score = 1
            if score == 0 and strand_pause == "-":
                score = -1

            if hashed not in wt_values.keys():
                if strand_pause == '+' and score > 0:
                    wt_values[hashed] = [abs(score)]
                elif strand_pause == '-' and score < 0:
                    wt_values[hashed] = [abs(score)]
                else:
                    wt_values[hashed] = [1]
            else:
                if strand_pause == '+' and score > 0:
                    wt_values[hashed].append(abs(score))
                elif strand_pause == '-' and score < 0:
                    wt_values[hashed].append(abs(score))
                else:
                    wt_values[hashed].append(1)

        for item in mutant_wig:
            score = item[key]

            if score == 0 and strand_pause == "+":
                score = 1
            if score == 0 and strand_pause == "-":
                score = -1

            if hashed not in mutant_values.keys():
                if strand_pause == '+' and score > 0:
                    mutant_values[hashed] = [abs(score)]
                elif strand_pause == '-' and score < 0:
                    mutant_values[hashed] = [abs(score)]
                else:
                    mutant_values[hashed] = [1]
            else:
                if strand_pause == '+' and score > 0:
                    mutant_values[hashed].append(abs(score))
                elif strand_pause == '-' and score < 0:
                    mutant_values[hashed].append(abs(score))
                else:
                    mutant_values[hashed].append(1)

    for key,value in pauses.items():
        hashed = str(key)+'_'+str(value)
        wt_scores = wt_values[hashed]
        mutant_scores = mutant_values[hashed]
        wt_mean = np.mean(wt_scores)
        mutant_mean = np.mean(mutant_scores)
        Log2 = round(np.log2(wt_mean/mutant_mean),2)
        pval = ttest_ind(wt_scores,mutant_scores, equal_var = False)[1]
        ranked_log2[hashed] = Log2
        ranked_pval[hashed] = pval
        raw_val[hashed] = [wt_mean,mutant_mean]

    return ranked_log2,ranked_pval,raw_val

def printer(outfile,pauses,genes,raw_vals,ranked_Log2,ranked_pval,ranked_FDR):
    with open(outfile,"w") as outp:
        outp.write('Coord,Strand,Gene,Dist,WT Mean,Mutant Mean,Log2FC,pval,FDR\n')

        for key,value in pauses.items():
            hashed = str(key)+'_'+str(value)
            pos = str(key)
            strand = str(value)
            wt_val = str(raw_vals[hashed][0])
            mutant_val =str(raw_vals[hashed][1])
            log2FC = str(ranked_Log2[hashed])
            pval = str(ranked_pval[hashed])
            qval = str(ranked_FDR[hashed])
            gene = genes[key][0]
            dist = genes[key][1]

            outp.write(pos+','+strand+','+gene+','+dist+','+wt_val+','+mutant_val+','+log2FC+','+pval+','+qval+'\n')


def main():
    parser = argparse.ArgumentParser(description='conducting differential pause stability analysis for construction of volcano plot')
    parser.add_argument('WT',type=str,help='Name of WT Condition Files')
    parser.add_argument('Mutant',type=str,help='Name of Mutant Conditon Files')
    parser.add_argument('outfile',type=str,help='name of output differential pause stability file')
    args = parser.parse_args()

    #load up wig files
    wt_wig = []
    for fyle in glob.glob(str(args.WT)+"*.wig"):
        wt_wig.append(read_wig(fyle))

    mutant_wig = []
    for fyle in glob.glob(str(args.Mutant)+"*.wig"):
        mutant_wig.append(read_wig(fyle))

    #build pause dict
    for fyle in glob.glob(str(args.WT)+"*.tab"):
        wt_pauses,wt_genes = read_tab(fyle)

    for fyle in glob.glob(str(args.Mutant)+"*.tab"):
        mutant_pauses,mutant_genes = read_tab(fyle)

    pauses = {**wt_pauses, **mutant_pauses}
    genes = {**wt_genes, **mutant_genes}

    #rank all pauses, obtain Log2FC and pval
    ranked_Log2,ranked_pval,raw_vals = pause_rank(pauses,wt_wig,mutant_wig)

    #clean up large dicts
    wt_wig.clear()
    mutant_wig.clear()

    #obtain FDR from pval
    sorted_ranked_pval = {k: v for k, v in sorted(ranked_pval.items(), key=lambda item: item[1])}
    pvals = [x[1] for x in sorted_ranked_pval.items()]
    FDR = multipletests(pvals, alpha=0.05, method='fdr_bh', is_sorted=True, returnsorted=False)[1]
    ranked_FDR = {}
    i = 0
    for item in sorted_ranked_pval:
        ranked_FDR[item] = FDR[i]
        i += 1

    #print everything out
    printer(args.outfile,pauses,genes,raw_vals,ranked_Log2,ranked_pval,ranked_FDR)

    #build graph
    data = pd.read_csv(args.outfile)
    visuz.gene_exp.volcano(d=data, lfc='Log2FC', pv='FDR',sign_line=True)

if __name__ == '__main__':
    main()
