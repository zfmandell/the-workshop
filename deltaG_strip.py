import sys
import glob

delta = []
for fyle in glob.glob('./CT/*.ct'):
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        try:
            delta.append(str(firstline.strip().split()[3]))
        except IndexError:
            delta.append('0')

with open(sys.argv[1],'w') as outp:
    outp.write('deltaG\n')
    for item in delta:
        outp.write(item+"\n")
