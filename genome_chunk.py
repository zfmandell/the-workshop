import sys
from Bio import SeqIO
import random

def read_fasta(genome_fasta):
    #reads in fasta
    fasta_sequences,fasta_dict =SeqIO.parse(open(genome_fasta),'fasta'),{}
    for fasta in fasta_sequences:
        fasta_dict['fasta'] = str(fasta.seq)
    return fasta_dict

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i + n]

def iter_sample_fast(iterable, samplesize):
    results = []
    iterator = iter(iterable)
    # Fill in the first samplesize elements:
    try:
        for _ in xrange(samplesize):
            results.append(iterator.next())
    except StopIteration:
        raise ValueError("Sample larger than population.")
    random.shuffle(results)  # Randomize their positions
    for i, v in enumerate(iterator, samplesize):
        r = random.randint(0, i)
        if r < samplesize:
            results[r] = v  # at a decreasing rate, replace random items
    return results



#chunked = [genome[i:i + 40] for i in xrange(0, len(genome), 40)]
#for i in xrange(0,1501):
    #return_dict[str(i)] =

with open(sys.argv[2],'w') as outp:
    i = 0
    for item in iter_sample_fast(chunks(read_fasta(sys.argv[1])['fasta'],40),1501):
        outp.write('>'+str(i)+"\n"+item+"\n")
        i+=1
