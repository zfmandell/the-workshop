import sys

final = []
with open(sys.argv[1],'r') as inp:
    for line in inp:
        end = int(line.strip().split('\t')[1])
        delta = float(line.strip().split('\t')[3])
        if delta > 0:
            strand = '+'
        else:
            strand = '-'
        final.append([end,strand,delta])

with open(sys.argv[2],'w') as outp:
    outp.write('end,strand,delta\n')
    for item in final:
        outp.write(str(item[0])+','+str(item[1])+','+str(item[2])+'\n')
