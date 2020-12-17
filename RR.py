from __future__ import division
import argparse

def read_gff3(gff3):
    #reads in fasta, returns dict of sequences, dict of strand info for gene
    return_dict,strand_dict = {},{}
    with open(gff3,'r') as inp:
        for line in inp:
            if line[0] != '#' and len(line.strip().split()) != 0:
                if str(line.strip().split()[2]) == 'gene':
                    return_dict[(line.strip().split()[8].split(";")[2].split("=")[1]).replace("-",".")] = [int(line.strip().split()[3]),int(line.strip().split()[4])]
                    strand_dict[(line.strip().split()[8].split(";")[2].split("=")[1]).replace("-",".")] = line.strip().split()[6]
    return return_dict,strand_dict

def read_UTR(UTR):
    return_dict = {}
    with open(UTR,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            return_dict[line.strip().split(",")[0]] = [int(line.strip().split(",")[1]),int(line.strip().split(",")[2])]
    return return_dict

def gff3_stitch(UTR,gff3,strands):
    return_dict = {}
    for key,value in UTR.iteritems():
        if strands[key] == '+':
            return_dict[key] = [gff3[key][0]-value[0],gff3[key][1]+value[1]]
        else:
            return_dict[key] = [gff3[key][1]+value[0],gff3[key][0]-value[1]]
    for key,value in gff3.iteritems():
        if key not in return_dict.keys():
            return_dict[key] = value
    return return_dict

def peak_reader(peaks):
    return_dict = {}
    with open(peaks,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            if line.strip().split(",")[1].split(":")[1] not in return_dict.keys():
                to_add = {}
                if line.strip().split(",")[2] != '-':
                    to_add['WT'] = [int(line.strip().split(",")[0])]
                else:
                    to_add['WT'] = []
                if line.strip().split(",")[3] != '-' and line.strip().split(",")[2] == '-':
                    to_add['KO'] = [int(line.strip().split(",")[0])]
                else:
                    to_add['KO'] = []
                return_dict[line.strip().split(",")[1].split(":")[1]] = to_add
            else:
                if line.strip().split(",")[2] != '-':
                    return_dict[line.strip().split(",")[1].split(":")[1]]['WT'].append(int(line.strip().split(",")[0]))
                if line.strip().split(",")[3] != '-' and line.strip().split(",")[2] == '-':
                    return_dict[line.strip().split(",")[1].split(":")[1]]['KO'].append(int(line.strip().split(",")[0]))
    return return_dict

def RR_reader(RR):
    with open(RR,'r') as inp:
        firstline = inp.readline()
        return_list = [line.strip().split(",")[0] for line in inp]
    return return_list

def KO_analysis(peaks,RR,genes,strands):
    KO,WT = [],[]
    for gene in RR:
        if gene in peaks.keys():
            if strands[gene] == '+':
                for peak in peaks[gene]['WT']:
                    WT.append(float(peak-genes[gene][0])/float(genes[gene][1]-genes[gene][0]))
                for peak in peaks[gene]['KO']:
                    KO.append(float(peak-genes[gene][0])/float(genes[gene][1]-genes[gene][0]))
            else:
                for peak in peaks[gene]['WT']:
                    WT.append(float(peak-genes[gene][0])/float(genes[gene][1]-genes[gene][0]))
                for peak in peaks[gene]['KO']:
                    KO.append(float(genes[gene][1]-peak)/float(genes[gene][1]-genes[gene][0]))
    return KO,WT

def list_bulk(a,b):
    if len(a) > len(b):
        while len(a) > len(b):
            b.append('NA')
    if len(b) > len(a):
        while len(b) > len(a):
            a.append('NA')
    return a,b

def main():
    parser = argparse.ArgumentParser(description='analysis of ratio of ratios for PNPase KO')
    parser.add_argument('RR',type=str,help='RR.csv')
    parser.add_argument('peaks',type=str,default='plus',help='peak csv file')
    parser.add_argument('UTR',type=str,default=None,help='<.txt> gene feature information')
    parser.add_argument('genbank',type=str,default=None,help='<.gff3> genbank file for BSub')
    parser.add_argument('outfile',type=str,default=None,help='name of output files')
    args = parser.parse_args()

    peaks = peak_reader(args.peaks)
    RR = RR_reader(args.RR)
    genes,strands = read_gff3(args.genbank)
    UTR = read_UTR(args.UTR)
    stitched = gff3_stitch(UTR,genes,strands)
    KO,WT = KO_analysis(peaks,RR,stitched,strands)
    KO_len,WT_len = list_bulk(KO,WT)
    KO_final = [x*100 for x in KO_len]
    WT_final = [x*100 for x in WT_len]

    with open(args.outfile,'w') as outp:
        outp.write('WT,KO\n')
        i = 0
        while i < len(WT_final):
            outp.write(str(WT_final[i])+","+str(KO_final[i])+"\n")
            i+=1

if __name__ == '__main__':
    main()
