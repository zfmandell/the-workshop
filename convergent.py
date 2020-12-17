import sys

with open(sys.argv[1],'r') as inp:
    lines = [x for x in inp.readlines() if x[0] != '#']
    for item in lines[:len(lines)-1]:
        if item.strip().split()[2] == 'CDS' and item.strip().split()[6] == '+' and lines[lines.index(item)+1].strip().split()[6] == '-':
            if 0 < int(lines[lines.index(item)+1].strip().split()[3]) - int(item.strip().split()[4]) <= 100:
                print '%d %d' % (int(item.strip().split()[4])-25,int(lines[lines.index(item)+1].strip().split()[3])+25)
            elif int(lines[lines.index(item)+1].strip().split()[3]) - int(item.strip().split()[4]) <= 0:
                print '%d %d' % (int(lines[lines.index(item)+1].strip().split()[3])-25,int(item.strip().split()[4])+25)
