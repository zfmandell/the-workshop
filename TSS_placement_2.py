import sys

GTF,TSS = [],{}

def opposite(stri):
    if stri == '+':
        return '-'
    else:
        return '+'

with open(sys.argv[1],'r') as inp:
    for line in inp:
        if line[0] != '@' and line.strip().split()[5] == '100M':
            if line.strip().split()[1] == '0':
                TSS[line.strip().split()[0]] = [int(line.strip().split()[3])+50,'+']
            elif line.strip().split()[1] == '16':
                TSS[line.strip().split()[0]] = [int(line.strip().split()[3])+50,'-']

with open(sys.argv[2],'r') as inp:
    for line in inp:
        if line.strip().split('\t')[3] in TSS.keys():
            if str(TSS[line.strip().split('\t')[3]][1]) == '+':
                ID = 'pacbio'
                src = 'GAU_Zach'
                feat = 'TSS'
                start = TSS[line.strip().split('\t')[3]][0]-1
                end = TSS[line.strip().split('\t')[3]][0]-1
                score = '.'
                strand = line.strip().split('\t')[6]
                phase = '.'
                attribute = line.strip().split('\t')[8]
                GTF.append([ID,src,feat,start,end,score,strand,phase,attribute])
            else:
                ID = 'pacbio'
                src = 'GAU_Zach'
                feat = 'TSS'
                start = TSS[line.strip().split('\t')[3]][0]-1
                end = TSS[line.strip().split('\t')[3]][0]-1
                score = '.'
                strand = opposite(line.strip().split('\t')[6])
                phase = '.'
                attribute = line.strip().split('\t')[8]
                GTF.append([ID,src,feat,start,end,score,strand,phase,attribute])

GTF_sorted = sorted(GTF, key=lambda x: x[3])

with open(sys.argv[3],'w') as outp:
    for item in GTF_sorted:
        outp.write('\t'.join(map(str,item)))
        outp.write('\n')
