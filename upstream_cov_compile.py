import glob
import sys

def translator(lyst,base):
    for item in base:
        if item[0]-3 <= lyst[0] <= item[0]+3 and lyst[1] == item[1][0]:
            return item[0]
            break
    return None

ends = {}
with open(sys.argv[1],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        ends[int(line.strip().split(",")[0])] = line.strip().split(",")[1]

final = {}
for fyle in sorted(glob.glob('*.csv')):
    strain = fyle.strip().split("_")[0]
    to_add = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            if line.strip().split(",")[1] == 'fwd':
                to_add[translator([int(line.strip().split(",")[0]),'+'],ends.items())] = ['+',line.strip().split(",")[4]]
            else:
                to_add[translator([int(line.strip().split(",")[0]),'-'],ends.items())] = ['-',line.strip().split(",")[4]]

    final[strain] = to_add

ends = {}
with open(sys.argv[1],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        ends[int(line.strip().split(",")[0])] = line.strip().split(",")[1]

with open(sys.argv[2],'w') as outp:
    outp.write('POT,WT,dA,dG,dR,dAG,dAR,dGR,dAGR\n')
    for key,value in ends.iteritems():
        outp.write(str(key)+",")
        if key in final['WT'].keys():
            if value == final['WT'][key][0]:
                outp.write(final['WT'][key][1]+",")
            else:
                outp.write('0,')
        else:
            outp.write('0,')

        if key in final['dA'].keys():
            if value == final['dA'][key][0]:
                outp.write(final['dA'][key][1]+",")
            else:
                outp.write('0,')
        else:
            outp.write('0,')

        if key in final['dG'].keys():
            if value == final['dG'][key][0]:
                outp.write(final['dG'][key][1]+",")
            else:
                outp.write('0,')
        else:
            outp.write('0,')

        if key in final['dR'].keys():
            if value == final['dR'][key][0]:
                outp.write(final['dR'][key][1]+",")
            else:
                outp.write('0,')
        else:
            outp.write('0,')

        if key in final['dAG'].keys():
            if value == final['dAG'][key][0]:
                outp.write(final['dAG'][key][1]+",")
            else:
                outp.write('0,')
        else:
            outp.write('0,')

        if key in final['dAR'].keys():
            if value == final['dAR'][key][0]:
                outp.write(final['dAR'][key][1]+",")
            else:
                outp.write('0,')
        else:
            outp.write('0,')

        if key in final['dGR'].keys():
            if value == final['dGR'][key][0]:
                outp.write(final['dGR'][key][1]+",")
            else:
                outp.write('0,')
        else:
            outp.write('0,')

        if key in final['dAGR'].keys():
            if value == final['dAGR'][key][0]:
                outp.write(final['dAGR'][key][1]+"\n")
            else:
                outp.write('0\n')
        else:
            outp.write('0\n')
