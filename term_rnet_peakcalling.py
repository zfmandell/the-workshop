#!/usr/bin/env python3

import argparse
import statistics
import numpy as np 

def term_reader(csv):
    return_dict = {}
    with open(csv,'r') as inp:
        firstline = inp.readline()
        for line in inp: 
            coord = int(line.strip().split(',')[0])
            strand = line.strip().split(',')[1]
            return_dict[coord] = strand
    
    return return_dict

def rnet_reader(bw):
    return_plus,return_minus = [],[]
    with open(bw,'r') as inp:
        for line in inp:
            count = int(line.strip().split('\t')[1])
            if count > 0:
                return_plus.append(count)
                return_minus.append(0)
            if count < 0:
                return_plus.append(0)
                return_minus.append(abs(count))
            elif count == 0:
                return_plus.append(0)
                return_minus.append(0)

    to_add = [0] * 101
    return_plus = to_add + return_plus + to_add
    return_minus = to_add + return_minus + to_add

    return return_plus,return_minus

def end_place(released,nascent,strand):
    ends,return_dict = [x for x,y in released.items() if y == strand],{}
    if strand == '+':
        for coord in ends:
            nascent_terminal = nascent[coord:coord+5]
            nascent_baseline = nascent[coord-150:coord+5]
            
            if all(v == 0 for v in nascent_terminal) or np.percentile(nascent_baseline,75) == 0:
                pass
            else:
                temp = [x/np.percentile(nascent_baseline,75) for x in nascent_terminal]
                return_dict[str(coord)+'_+'] = temp
                
    else:
        for coord in ends:
            nascent_terminal = nascent[coord-4:coord+1][::-1]
            nascent_baseline = nascent[coord-4:coord+151]

            if all(v == 0 for v in nascent_terminal) or np.percentile(nascent_baseline,75) == 0:
                pass
            else:
                temp = [x/np.percentile(nascent_baseline,75) for x in nascent_terminal]
                return_dict[str(coord)+'_-'] = temp
    
    return return_dict

def end_cruncher(ends,percentile):
    threshold = 2 ** percentile
    return_list = []
    for key,value in ends.items():
        if max(value) >= percentile:
            end_index = value.index(max(value))
            return_list.append([int(key.split('_')[0]),key.split('_')[1],end_index])

    return return_list


def main():
    parser = argparse.ArgumentParser(description='determines 3-OH end idendity by combining Term-seq and RNET-seq information')
    parser.add_argument('Termseq',type=str,help='list of Term-seq called 3-OH ends: <.csv>, no default value, must be inputted')
    parser.add_argument('RNETseq',type=str,help='list of RNET-seq called 3-OH ends: <.bw>, no default value, must be inputted')

    parser.add_argument('output',type=str,help='name of output file')
   
    args = parser.parse_args()

    released = term_reader(args.Termseq)
    nascent_plus,nascent_minus = rnet_reader(args.RNETseq)
    authentic_plus,authentic_minus = end_place(released,nascent_plus,'+'),end_place(released,nascent_minus,'-')
    authentic = {**authentic_plus,**authentic_minus}
    authentic_crunched = end_cruncher(authentic,np.percentile([np.log2(max(item)) for item in authentic.values()], 25))
    authentic_ends = sorted(authentic_crunched)
    

    with open(args.output,'w') as outp:
        outp.write('coord,strand,dist\n')
        for item in authentic_ends:
            outp.write(str(item[0])+','+str(item[1])+','+str(item[2])+'\n')

if __name__ == '__main__':
    main()
