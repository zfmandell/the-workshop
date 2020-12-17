import sys

with open(sys.argv[1],'r') as inp:
    with open('inverse_'+sys.argv[1],'w') as outp:
        for line in inp:
            outp.write(line.strip().split('\t')[0]+'\t'+line.strip().split('\t')[1]+'\t'+line.strip().split('\t')[2]+'\t'+'-'+line.strip().split('\t')[3]+'\n')
