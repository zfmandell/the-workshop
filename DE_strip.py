import sys

with open(sys.argv[1],'r') as inp:
    i,q = 0,0
    firstline = inp.readline()
    for line in inp:
        try:
            if float(line.strip().split(',')[2]) >= 2 and float(line.strip().split(',')[6]) <= 0.005:
                i+=1
            elif float(line.strip().split(',')[2]) <= -2 and float(line.strip().split(',')[6]) <= 0.005:
                q+=1
        except ValueError:
            pass

print 'lower'
print q
print 'higher'
print i
