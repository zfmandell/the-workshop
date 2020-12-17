from __future__ import division
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from itertools import islice
import collections
from natsort import natsorted

def window(seq, n=2):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result

YY = {}
for record in SeqIO.parse(sys.argv[1], "fasta"):
    q=0
    for item in window(str(record.seq)):
        if ''.join(item) == 'UC' or ''.join(item) == 'CC':
            if str(q) in YY.keys():
                YY[str(q)] +=1
                q+=1
            else:
                YY[str(q)] = 1
                q+=1
        else:
            q+=1

num = len([1 for line in open(sys.argv[1]) if line.startswith(">")])

with open(sys.argv[2],'w') as outp:
    outp.write("pos,#\n")
    for key in natsorted(YY):
        outp.write(str(key)+","+str(YY[key]/num)+"\n")
