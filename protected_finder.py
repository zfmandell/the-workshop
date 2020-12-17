#!/usr/bin/env python3

import argparse
import numpy as np
from scipy.signal import find_peaks
from scipy.signal import argrelextrema

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

def read_coverage(cov_fyle):
    #create list of coverage values
    return_list = []
    with open(cov_fyle,"r") as inp:
        for line in inp:
            return_list.append(int(line.strip().split("\t")[2]))
    return return_list

def contiguous_regions(condition):
    """Finds contiguous True regions of the boolean array "condition". Returns
    a 2D array where the first column is the start index of the region and the
    second column is the end index."""

    # Find the indicies of changes in "condition"
    d = np.diff(condition)
    idx, = d.nonzero()

    # We need to start things after the change in "condition". Therefore,
    # we'll shift the index by 1 to the right.
    idx += 1

    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = np.r_[0, idx]

    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = np.r_[idx, condition.size] # Edit

    # Reshape the result into two columns
    idx.shape = (-1,2)
    return idx

def main():
    parser = argparse.ArgumentParser(description='determines positions and lengths of DNAseq coverage')
    parser.add_argument('DNAseq',type=str,help='file name of DNAseq cov file')
    parser.add_argument('RNET_fwd',type=str,help='file name of RNETseq fwd cov file')
    parser.add_argument('RNET_rev',type=str,help='file name of RNETseq rev cov file')
    parser.add_argument('fasta',type=str,help='fasta file of genome')
    parser.add_argument('output',type=str,help='output base name')
    parser.add_argument('-threshold_ratio',type=int,default=10,help='threshold of DNAseq:RNETseq coverage to be considered differential, default = 10')
    parser.add_argument('-threshold_percentile',type=int,default=25,help='percentile of DNAseq coverage to be considered a peak, default = 25')
    parser.add_argument('-threshold_length',type=int,default=25,help='minimum length of protected region, default = 25')


    args = parser.parse_args()

    seq = genome_yield(args.fasta)

    cov_DNA = np.array(read_coverage(args.DNAseq))
    cov_RNET_fwd = np.array(read_coverage(args.RNET_fwd))
    cov_RNET_rev = np.array(read_coverage(args.RNET_rev))
    condition = np.abs(cov_DNA) > np.percentile(cov_DNA,args.threshold_percentile)


    return_peaks = []
    for start, stop in contiguous_regions(condition):
        segment_DNA = cov_DNA[start:stop]
        segment_RNET_fwd = cov_RNET_fwd[start:stop]
        segment_RNET_rev = cov_RNET_rev[start:stop]
        segment_RNET = [a + b for a,b in zip(list(segment_RNET_fwd),list(segment_RNET_rev))]

        DNA = np.sum(segment_DNA)/(stop-start)
        RNET = sum(segment_RNET)/(stop-start)

        if DNA > RNET*args.threshold_ratio and stop-start > args.threshold_length:
            if np.sum(segment_RNET_fwd) > np.sum(segment_RNET_rev):
                return_peaks.append([start,stop,'+',stop-start,DNA,RNET,DNA/RNET])
            else:
                return_peaks.append([start,stop,'-',stop-start,DNA,RNET,DNA/RNET])

    with open(args.output+'.csv','w') as outp:
        outp.write('start,stop,strand,length,DNA,RNET,DNA:RNET\n')
        for item in return_peaks:
            outp.write(str(item[0])+","+str(item[1])+","+str(item[2])+","+str(item[3])+","+str(item[4])+","+str(item[5])+","+str(item[6])+'\n')

    with open(args.output+'.fa','w') as outp:
        for item in return_peaks:
            if item[2] == '+':
                outp.write('>'+str(str(item[0])+"_"+str(item[1]))+'\n')
                outp.write(seq[item[0]:item[1]+1]+'\n')
            else:
                outp.write('>'+str(str(item[0])+"_"+str(item[1]))+'\n')
                outp.write(reverse_complement(seq[item[0]:item[1]+1])+'\n')


if __name__ == '__main__':
    main()
