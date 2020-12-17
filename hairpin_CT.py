#!/usr/bin/env python

import sys
import glob

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i

def collect_information(directory):
    #Snags a directory worth of CT files,runs terminal report on each CT file, outputs list of terminal reports
    information = []
    for fyle in glob.glob(directory+'/*.ct'):
        information.append(terminal_report(fyle))
    return information

def terminal_report(fyle):
    #generates comprehensive information list of lists, indexed as mentioned above (#1-5)
    return_list = []
    if file_len(fyle) == 41:
        with open(fyle,'r') as inp:
            firstline = inp.readline()
            if len(firstline.strip().split()) > 2:
                delta,name,temp,index = firstline.strip().split()[3],firstline.strip().split()[4],[],0
                for line in inp:
                    temp.append(int(line.strip().split()[4]))
                for item in reversed(temp):
                    if item == 0:
                        index+=1
                    else:
                        break
                return_list.append([name.split('-')[0],name.split('-')[1],name.split('-')[2],delta,index])
            else:
                return_list.append([firstline.strip().split()[1].split('-')[0],firstline.strip().split()[1].split('-')[1],firstline.strip().split()[1].split('-')[2],0,'NA'])
    else:
        with open(fyle,'r') as inp:
            fyle = inp.readlines()
            chunked,to_probe = [fyle[i:i + 42] for i in xrange(0, len(fyle), 42)],[]
            for item in chunked:
                temp_list,temp_index = [],0
                for sub_item in item[1:]:
                    temp_list.append(int(sub_item.strip().split()[4]))
                for sub_item in reversed(temp_list):
                    if sub_item == 0:
                        temp_index+=1
                    else:
                        break
                to_probe.append([item[0].strip().split()[4],float(item[0].strip().split()[3]),temp_index])

            lowest = min(x[2] for x in to_probe)
            finalists = []
            for item in to_probe:
                if item[2] == lowest:
                    finalists.append(item)
            if len(finalists) > 1:
                final = []
                highest = min(x[1] for x in finalists)
                for item in finalists:
                    if item[1] == highest:
                        final.append(item)
                flattened = [item for sublist in final for item in sublist]
                return_list.append([flattened[0].split("-")[0],flattened[0].split("-")[1],flattened[0].split("-")[2],flattened[1],flattened[2]])
            else:
                flattened = [item for sublist in finalists for item in sublist]
                return_list.append([flattened[0].split("-")[0],flattened[0].split("-")[1],flattened[0].split("-")[2],flattened[1],flattened[2]])

    return return_list

directory = sys.argv[1]

with open(sys.argv[2],'w') as outp:
    outp.write('gene,coordinate,score,deltaG,pos\n')
    for item in collect_information(directory):
        flattened = [item for sublist in item for item in sublist]
        if flattened[4] == 'NA':
            pass
        elif int(flattened[4])+10 > 12:
            pass
        elif float(flattened[3]) > -3.0:
            pass
        else:
            outp.write(str(flattened[0])+","+str(flattened[1])+","+str(flattened[2])+","+str(flattened[3])+","+str(int(flattened[4])+10)+"\n")
