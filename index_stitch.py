#!/usr/bin/env python
from __future__ import division
import argparse
import glob
import operator

def collect_CT_info(CT_PATH):
    #Snags a directory worth of CT files,runs terminal report on each CT file, outputs list of terminal reports

    information,return_dict = [],{}
    for fyle in glob.glob(str(CT_PATH)+'*.ct'):
        information.append(terminal_report(fyle))

    ends = [float(item[0]) for item in information]
    start = [item[1] for item in information]

    q=0
    for item in ends:
        return_dict[item] = start[q]
        q+=1
    return return_dict

def terminal_report(fyle):
    #generates comprehensive information list of lists, indexed as mentioned above (#1-5)
    struct,return_dict = [],{}

    with open(fyle,'r') as inp:
        firstline = inp.readline()
        try:
            end = firstline.strip().split()[4].split("_")[0].split(":")[1]
        except IndexError:
            end = firstline.strip().split()[1].split("_")[0].split(":")[1]

        struct = [int(x.split()[4]) for x in [next(inp) for x in xrange(61)]]
        rev_struct = list(reversed(struct))

        q = 0
        for item in rev_struct:
            q+=1
            if int(item) > 0:
                distance_start = q
                break

        try:
            return [end,distance_start]
        except UnboundLocalError:
            return [end,'NA']

def read_gff3(gff3):
    #reads in fasta, returns dict of sequences, dict of strand info for gene
    return_dict,strand_dict = {},{}
    with open(gff3,'r') as inp:
        for line in inp:
            if line[0] != '#' and len(line.strip().split()) != 0:
                if str(line.strip().split()[2]) == 'gene':
                    return_dict[(line.strip().split()[8].split(";")[2].split("=")[1]).replace("-",".")] = [int(line.strip().split()[3]),int(line.strip().split()[4])]
                    strand_dict[(line.strip().split()[8].split(";")[2].split("=")[1]).replace("-",".")] = line.strip().split()[6]
    return sorted(return_dict.items(), key=operator.itemgetter(1)),strand_dict,return_dict

def read_MFE(MFE):
    return_dict = {}
    with open(MFE,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            return_dict[float(line.strip().split(",")[0].split(":")[1])] = float(line.strip().split(",")[1])
    return return_dict

def read_ends(ends):
    return_list = []
    with open(ends,'r') as inp:
        for line in inp:
            return_list.append(float(line.strip()))
    return return_list

def read_UTR(UTR):
    return_dict = {}
    with open(UTR,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            return_dict[line.strip().split(",")[0]] = [int(line.strip().split(",")[1]),int(line.strip().split(",")[2])]
    return return_dict

def read_bedgraph(bedgraph):
    strand_dict,delta_dict = {},{}
    with open(bedgraph,'r') as inp:
        for line in inp:
            if float(line.strip().split("\t")[3]) < 0:
                strand_dict[float(line.strip().split("\t")[2])] = '-'
            elif float(line.strip().split("\t")[3]) > 0:
                strand_dict[float(line.strip().split("\t")[2])] = '+'
            delta_dict[float(line.strip().split("\t")[2])] = abs(float(line.strip().split("\t")[3]))
    return strand_dict,delta_dict

def gene_compile(ends,seqs,gene_strands,UTR,end_strands,raw_dict):
    #returns dict where coordinate : genomical_feature (either CDS, or between two CDS)
    feature_dict,return_dict = {},{}
    for item in ends:
        for seq in seqs:
            if item < seqs[0][1][0]:
                feature_dict[item] = "-".join([seqs[0][0],seqs[-1][0]])
            elif seq[1][0] <= item <= seq[1][1]:
                if gene_strands[seq[0]] == end_strands[item]:
                    feature_dict[item] = 'sense_CDS:'+str(seq[0])
                else:
                    feature_dict[item] = 'antisense_CDS:'+str(seq[0])
            elif seqs.index(seq) < len(seqs) - 1:
                if seq[1][1] < item < seqs[seqs.index(seq)+1][1][0]:
                    feature_dict[item] = "-".join([seq[0],seqs[seqs.index(seq)+1][0]])
            elif seqs.index(seq) == len(seqs) - 1:
                if seq[1][1] < item:
                    feature_dict[item] = "-".join([seq[0],seqs[0][0]])

    for key,value in feature_dict.iteritems():
        if len(value.split("-")) == 2:
            if gene_strands[value.split("-")[0]] and gene_strands[value.split("-")[1]] == '+':
                if end_strands[key] == '-':
                    return_dict[key] ='antisense:'+str(value)
                else:
                    if raw_dict[value.split("-")[0]][1] < raw_dict[value.split("-")[1]][0]:
                        left = value.split("-")[0]
                        right = value.split("-")[1]
                    else:
                        left = value.split("-")[1]
                        right = value.split("-")[0]

                    if key < raw_dict[left][1]+UTR[left][1] and key > raw_dict[right][0]-UTR[right][0]:
                        return_dict[key] = 'FPUTR:'+str(right)+"_"+"TPUTR:"+str(left)
                    elif key < raw_dict[left][1]+UTR[left][1]:
                        return_dict[key] = "TPUTR:"+str(left)
                    elif  key > raw_dict[right][0]-UTR[right][0]:
                        return_dict[key] = 'FPUTR:'+str(right)
                    else:
                        return_dict[key] = value

            elif gene_strands[value.split("-")[0]] and gene_strands[value.split("-")[1]] == '-':
                if end_strands[key] == '+':
                    return_dict[key] ='antisense:'+str(value)
                else:

                    if raw_dict[value.split("-")[0]][0] > raw_dict[value.split("-")[1]][1]:
                        left = value.split("-")[0]
                        right = value.split("-")[1]
                    else:
                        left = value.split("-")[1]
                        right = value.split("-")[0]

                    if key > raw_dict[left][0]-UTR[left][1] and key < raw_dict[right][1]+UTR[right][0]:
                        return_dict[key] = 'FPUTR:'+str(right)+"_"+"TPUTR:"+str(left)
                    elif key > raw_dict[left][0]-UTR[left][1]:
                        return_dict[key] = "TPUTR:"+str(left)
                    elif  key < raw_dict[right][1]+UTR[right][0]:
                        return_dict[key] = 'FPUTR:'+str(right)
                    else:
                        return_dict[key] = value

            elif gene_strands[value.split("-")[0]] != gene_strands[value.split("-")[1]]:
                if end_strands[key] == gene_strands[value.split("-")[0]]:
                    if end_strands[key] == "+":
                        if key > raw_dict[value.split("-")[0]][0]-UTR[value.split("-")[0]][0]:
                            return_dict[key] = 'FPUTR:'+str(value.split("-")[0])
                        elif key < raw_dict[value.split("-")[0]][0]+UTR[value.split("-")[0]][1]:
                            return_dict[key] = 'TPUTR:'+str(value.split("-")[0])
                        else:
                            return_dict[key] = value
                    else:
                        if key < raw_dict[value.split("-")[0]][1]+UTR[value.split("-")[0]][0]:
                            return_dict[key] = 'FPUTR:'+str(value.split("-")[0])
                        elif key < raw_dict[value.split("-")[0]][0]-UTR[value.split("-")[0]][1]:
                            return_dict[key] = 'TPUTR:'+str(value.split("-")[0])
                        else:
                            return_dict[key] = value
                else:
                    if end_strands[key] == "+":
                        if key > raw_dict[value.split("-")[1]][0]-UTR[value.split("-")[1]][0]:
                            return_dict[key] = 'FPUTR:'+str(value.split("-")[1])
                        elif key < raw_dict[value.split("-")[1]][0]+UTR[value.split("-")[1]][1]:
                            return_dict[key] = 'TPUTR:'+str(value.split("-")[1])
                        else:
                            return_dict[key] = value
                    else:
                        if key < raw_dict[value.split("-")[1]][1]+UTR[value.split("-")[1]][0]:
                            return_dict[key] = 'FPUTR:'+str(value.split("-")[1])
                        elif key < raw_dict[value.split("-")[1]][0]-UTR[value.split("-")[1]][1]:
                            return_dict[key] = 'TPUTR:'+str(value.split("-")[1])
                        else:
                            return_dict[key] = value

        else:
            return_dict[key] = value
    return sorted(return_dict.items(), key=operator.itemgetter(0))

def writer(compiled,MFE,delta,CT,outfile):
    with open(outfile,'w') as outp:
        outp.write("end_coordinate,feature,deltaG,delta,struct_dist\n")
        for item in compiled:
            outp.write(str(item[0])+","+str(item[1])+","+str(MFE[item[0]])+","+str(delta[item[0]])+","+str(CT[item[0]])+"\n")

def main():
    parser = argparse.ArgumentParser(description='will create .csv file of all 3 prime ends and information for each')
    parser.add_argument('ends',type=str,default=None,help='<.csv> list of 3 prime ends')
    parser.add_argument('MFE',type=str,default=None,help='<.csv> list of MFE values')
    parser.add_argument('UTR',type=str,default=None,help='<.txt> gene feature information')
    parser.add_argument('genbank',type=str,default=None,help='<.gff3> genbank file for BSub')
    parser.add_argument('bedgraph',type=str,default=None,help='<.bedgraph> bedgraph file of peaks')
    parser.add_argument('CT_path',type=str,default=None,help='path to CT directory')
    parser.add_argument('outfile',type=str,default=None,help='name of output file')
    args = parser.parse_args()

    seqs,gene_strands,raw_dict = read_gff3(args.genbank)
    end_strands,end_delta = read_bedgraph(args.bedgraph)
    MFE,ends,UTR,CT = read_MFE(args.MFE),read_ends(args.ends),read_UTR(args.UTR),collect_CT_info(args.CT_path)
    compiled = gene_compile(ends,seqs,gene_strands,UTR,end_strands,raw_dict)

    writer(compiled,MFE,end_delta,CT,args.outfile)

if __name__ == '__main__':
    main()
