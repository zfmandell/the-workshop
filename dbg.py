import collections, sys
from Bio import Seq, SeqIO, SeqRecord

def twin(km):
    return Seq.reverse_complement(km)

def kmers(seq,k):
    for i in xrange(len(seq)-k+1):
        yield seq[i:i+k]

def fw(km):
    for x in 'ACGT':
        yield km[1:]+x

def bw(km):
    for x in 'ACGT':
        yield x + km[:-1]

def build(fn,k=3,limit=0):
    d = collections.defaultdict(int)

    for f in fn:
        reads = SeqIO.parse(f,'fastq')
        for read in reads:
            seq_s = str(read.seq)
            seq_l = seq_s.split('N')
            for seq in seq_l:
                for km in kmers(seq,k):
                    d[km] +=1
                seq = twin(seq)
                for km in kmers(seq,k):
                    d[km] += 1

def get_contig_forward(d,km):
    c_fw = [km]

    while True:
        if sum(x in d for x in fw(c_fw[-1])) != 1:
            break

        cand = [x for x in fw(c_fw[-1]) if x in d][0]
        if cand in c_fw:
            break

        if sum(x in d for x in bw(cand)) != 1:
            break

        c_fw.append(cand)

    return c_fw

def get_contig(d,km):
    c_fw = get_contig_forward(d,km)
    c_bw = get_contig_forward(d,twin(km))

    c = [twin(x) for x in c_bw[-1:0:-1]] + c_fw
    s = c[0] + ''.join(x[-1] for x in c[1:])
    return s,c

def all_contigs(d):
    done = set()
    r = []
    for x in d:
        if x not in done:
            s,c = get_contig(d,x)
            for y in c:
                done.add(y)
                done.add(twin(y))
            r.append(s)
    return r

def print_dbg(cs):
    for i,x in enumerate(cs):
        print('>contig%d\n%s\n'%(i,x))

if __name__ == "__main__":
    k = int(sys.argv[1])
    d = build(sys.argv[2:],k,0)
    cs = all_contigs(d)
    print_dbg(cs)
