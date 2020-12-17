from __future__ import division
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from itertools import islice
import collections
from natsort import natsorted

def window(seq,n):
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
    for item in [str(record.seq)[i:i+20] for i in xrange(0, len(str(record.seq))-(20-1), 10)]:
        if 'UC' in item and 'CC' not in item:
            if str(q) in YY.keys():
                YY[str(q)] += item.count("UC")

            else:
                YY[str(q)] = item.count("UC")

            q+=1
        elif 'CC' in item and 'UC' not in item:
            if str(q) in YY.keys():
                YY[str(q)] += item.count("CC")

            else:
                YY[str(q)] = item.count("CC")
            q+=1
        elif 'CC' in item and 'UC' not in item:
            if str(q) in YY.keys():
                YY[str(q)] += (item.count("UC") + item.count("CC"))

            else:
                YY[str(q)] = (item.count("UC") + item.count("CC"))
            q+=1

num = len([1 for line in open(sys.argv[1]) if line.startswith(">")])

with open(sys.argv[2],'w') as outp:
    outp.write("window,#\n")
    for key in natsorted(YY):
        outp.write(str(key)+","+str(YY[key]/num)+"\n")
