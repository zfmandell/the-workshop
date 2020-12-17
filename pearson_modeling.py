from __future__ import division
import sys
from scipy.stats.stats import pearsonr
from itertools import islice
from Bio import SeqIO

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

def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def is_list(s):
    if isinstance(s, list) == True:
        return True
    else:
        return False

def remake_rtsc(rtsc_fyle,fasta_fyle):
    handle_rtsc = read_in_rtsc(rtsc_fyle)
    handle_fasta = fasta_fyle

    output_dict = {}

    for record in SeqIO.parse(handle_fasta, "fasta"):

        return_stops = []
        ident = str(record.id)
        starting_seq = str(record.seq)
        starting_stops = handle_rtsc[ident]

        i = 0

        while i < len(str(record.seq))-1:
            if str(starting_seq[i]).lower() == 'a' or str(starting_seq[i]).lower() == 'c':
                return_stops.append(starting_stops[i])

            i += 1

        output_dict[ident] = return_stops
    return output_dict

rtsc = [sys.argv[1],sys.argv[2]]
transcripts = sys.argv[3]
coverage = sys.argv[4]
temp = sys.argv[5]
fasta = sys.argv[6]
transcript_list = []

with open(transcripts,'r') as trans:
    for item in trans:
        item = item.strip()
        transcript_list.append(item)

with open(coverage,'r') as cov:
    coverage_dict = {}
    firstline = cov.readline()
    for item in cov:
        item = item.strip()
        contents = item.split(",")
        mean = (float(contents[1])+float(contents[2]))/2
        difference = abs(float(contents[1])-float(contents[2]))
        to_add = [mean,difference]
        coverage_dict[contents[0]] = to_add

correl = {}
spec_one,spec_two = remake_rtsc(rtsc[0],fasta),remake_rtsc(rtsc[1],fasta)

for key, value in spec_one.iteritems():
    if key in transcript_list:
        correl[key] = (pearsonr(value,spec_two[key]))


default_name = '_'.join(sorted([x.replace('.rtsc','') for x in [sys.argv[1],sys.argv[2]]]))+'_modeling_correlation.csv'

with open(default_name,"a+") as outp:
    outp.write("transcript,correlation,mean_coverage,coverage_difference,length,temp")
    outp.write("\n")
    for key, value in correl.iteritems():
        if value[0] > 0:
            outp.write(str(key)+","+str(value[0])+",")

            for key_cov, value_cov in coverage_dict.iteritems():

                if key_cov == key:
                    outp.write(str(value_cov[0])+","+str(value_cov[1])+",")

            length = 0
            contents = key.split("|")

            new_contents = []

            for item in contents:
                if ":" in str(item):
                    new_contents.append(item.split(":"))
                else:
                    new_contents.append(item)

            for item in new_contents:
                if is_list(item) == True:
                    for sub_item in item:
                        if is_int(sub_item) == True:
                            length += int(sub_item)
                else:
                    if is_int(item) == True:
                        length += int(item)

            outp.write(str(length)+",")
            outp.write(str(temp))
            outp.write("\n")
