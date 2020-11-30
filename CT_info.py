from __future__ import division
import glob
import sys


def collect_infomation():
    #Snags a directory worth of CT files,runs terminal report on each CT file, outputs list of terminal reports
    information = []
    for fyle in glob.glob('*.ct'):
        information.append(terminal_report(fyle))
    return information

def terminal_report(fyle):
    #generates comprehensive information list of lists, indexed as mentioned above (#1-5)
    struct = []

    with open(fyle,'r') as inp:
        firstline = inp.readline()
        try:
            end = firstline.strip().split()[4].split("_")[0].split(":")[1]
        except IndexError:
            end = firstline.strip().split()[1].split("_")[0].split(":")[1]

        struct = [int(x.split()[4]) for x in [next(inp) for x in xrange(61)]]

        rev_struct = list(reversed(struct))

        q = 0
        for item in rev_struct:
            q+=1
            if int(item) > 0:
                distance_start = q
                break

        struct_tab = 0
        for item in rev_struct:
            if item > 0:
                struct_tab+=1
        struct = (struct_tab/61)*100

        try:
            return [end,distance_start,struct]
        except UnboundLocalError:
            return [end,0,struct]

information = collect_infomation()
ends = [item[0] for item in information]
start = [item[1] for item in information]
struct = [item[2] for item in information]

with open(sys.argv[1],'w') as outp:
    q = 0
    while q < len(ends)-1:
        outp.write(str(ends[q])+","+str(start[q])+","+str(struct[q])+"\n")
        q+=1
