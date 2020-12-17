from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import translate
import commands
import sys
import subprocess
import shlex
import os
import re
import urllib2

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases

def read_fasta(genome_fasta):
    fasta_sequences,fasta_dict =SeqIO.parse(open(genome_fasta),'fasta'),{}
    for fasta in fasta_sequences:
        fasta_dict[fasta.id] = str(fasta.seq)
    return fasta_dict

def mod_string(s):
    final = s # we need to add the original string initially.
    i = 0
    while len(final) % 3 != 0:
        final += 'N' # keep adding char at alternate indexes
        if i == len(s)-1: # if we are at the end of the string, reset
            i = 0
        else:
            i += 1 # else move to next char
    return final

accessions = read_fasta(sys.argv[1])

to_remove = []
for key,value in accessions.iteritems():
    if key.split("|")[1] == 'MULTISPECIES:':
        to_remove.append(key)
for item in to_remove:
    del accessions[item]

Entrez.api_key = "insert NCBI API KEY HERE"
Entrez.email = "zxm44@psu.edu"
final,species_total = {},[]
for key,value in accessions.iteritems():
    if key.split("|")[-2] == 'sp.':
        species = "_".join(key.split("|")[-3:]).strip('[').strip(']')
    else:
        species = "_".join(key.split("|")[-2:]).strip('[').strip(']')
    if species not in species_total:
        with open('./temp_query.fa','w') as outp:
            outp.write('>'+key+"\n"+value)
        genome = commands.getstatusoutput('elink -db protein -id '+ str(key.split("|")[0])+' -target nuccore|efetch -format acc')[1]
        try:
            handle = Entrez.efetch(db="nucleotide", id=genome, rettype="fasta", retmode="text")
        except urllib2.HTTPError as err:
            if err.code == 400:
                continue
            else:
                raise
        with open('temp_db.fa','w') as outp:
            outp.write(handle.read())
        for record in SeqIO.parse('temp_db.fa', "fasta"):
            sequence = record.seq
        subprocess.call(shlex.split('makeblastdb -dbtype nucl -in ./temp_db.fa  -out ./temp_db'))
        result = commands.getstatusoutput('tblastn -db ./temp_db -query temp_query.fa -outfmt "6 pident sstart send sstrand"| head -n1')[1].split()
        os.system('rm temp*')
        if float(result[0]) == 100.000 and int(result[1]) > int(result[2]):
            battle = reverse_complement(sequence[int(result[1])+1:int(result[1])+302])
            for i in xrange(0,3):
                tester = translate(mod_string(battle[i:])).split('*')
                q=0
                for item in tester:
                    if re.search(r'[K|R].[K|R]$',item) != None:
                        if 'M' in item:
                            final[species] = item[item.rfind('M'):]
                        elif 'M' not in item and 'V' in item:
                            final[species] = item[item.rfind('V'):]
                        elif 'M' not in item and 'V' not in item and 'L' in item:
                            final[species] = item[item.rfind('L'):]
                q+=1
        elif float(result[0]) == 100.000 and int(result[2]) > int(result[1]):
            battle = str(sequence[int(result[1])-301:int(result[1])])
            for i in xrange(0,3):
                tester = translate(mod_string(battle[i:])).split('*')
                q=0
                for item in tester:
                    if re.search(r'[K|R].[K|R]$',item) != None:
                        if 'M' in item:
                            final[species] = item[item.rfind('M'):]
                        elif 'M' not in item and 'V' in item:
                            final[species] = item[item.rfind('V'):]
                        elif 'M' not in item and 'V' not in item and 'L' in item:
                            final[species] = item[item.rfind('L'):]
                q+=1
        species_total.append(species)

with open(sys.argv[2],'w') as outp:
    for key,value in final.iteritems():
        if len(value) >= 7 and 'sp.' not in key and hasNumbers(key) == False and len(value) <= 15::
            outp.write(">"+str(key)+"\n"+str(value)+"\n")
