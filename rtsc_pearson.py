from __future__ import division
import sys
import math
from itertools import islice
from scipy.stats.stats import pearsonr
import numpy as np
from Bio import SeqIO

def mean(numbers):
    return sum(numbers) / len(numbers)

def gen_dict(dict_1,dict_2):
    return_dict = {}
    for key, value in dict_1.iteritems():
        return_dict[key] = (pearsonr(value,dict_2[key]))
    return return_dict

def read_in_rtsc(rtsc_fyle):
    '''Reads a <.rtsc> file into a dictionary, transcript_name:[list of stop numbers]'''
    information = {}
    with open(rtsc_fyle,'r') as f:
        while True:
            next_n_lines = list(islice(f, 3))
            if not next_n_lines:
                break
            transcript,stops,empty_line = [n.strip() for n in next_n_lines]
            information[transcript] = [int(x) for x in stops.split('\t')]
    return information

def remake_rtsc(rtsc_fyle,fasta_fyle):
    handle_rtsc = rtsc_fyle
    handle_fasta = fasta_fyle

    output_dict = {}

    for record in SeqIO.parse(handle_fasta, "fasta"):

        return_stops = []
        ident = str(record.id)
        starting_seq = str(record.seq)
        starting_stops = handle_rtsc[ident]

        i = 0

        while i < len(str(record.seq))-1:
            if str(starting_seq[i]).lower() == 'g':
                return_stops.append(starting_stops[i])

            i += 1

        output_dict[ident] = return_stops
    return output_dict

def pearson(replicate_1,replicate_2):
    #pass empty dict and array (cor)
    temp_array = []
    mean_array = []
    for key, value in replicate_1.iteritems():
        temp_array.append(spearmanr(value,replicate_2[key]))

    for item in temp_array:
        if isinstance(item[0], float) == True:
            mean_array.append(item[0])
    return np.nanmean(mean_array)

rtsc = [sys.argv[1],sys.argv[2]]
transcripts = sys.argv[3]
fasta= sys.argv[4]
transcript_list = []

with open(transcripts,'r') as trans:
    for item in trans:
        item = item.strip()
        transcript_list.append(item)

correl = {}
replicate_1,replicate_2 = read_in_rtsc(rtsc[0]),read_in_rtsc(rtsc[1])
mod_1,mod_2 = remake_rtsc(replicate_1,fasta),remake_rtsc(replicate_2,fasta)

for key, value in mod_1.iteritems():
    if key in transcript_list:
        correl[key] = (pearsonr(value,mod_2[key]))



default_name = '_'.join(sorted([x.replace('.rtsc','') for x in [sys.argv[1],sys.argv[2]]]))+'_pearson_correlation.csv'
with open(default_name,"a+") as outp:
    outp.write("transcript,pearson_correlation,pvalue")
    outp.write("\n")
    for key, value in correl.iteritems():
        if value[0] > 0:

            outp.write(str(key)+","+str(value[0])+","+str(value[1]))
            outp.write("\n")
