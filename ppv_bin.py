#!/usr/bin/env python

import os,sys
from collections import defaultdict
import subprocess
import shutil
import distutils.dir_util
import argparse


def is_float(item):
    try:
        float(value)
    except ValueError:
        return False
    except TypeError:
        return False

def main():
    parser = argparse.ArgumentParser(description='binning and grabbing ps of interest, MUST RUN IN FOLDER WITH PS AND CT')
    parser.add_argument('ppv',default = None, help = '<.txt > ppv file')
    parser.add_argument('condition',default=None,help = '<.txt> biological condition')

    args = parser.parse_args()

    with open(args.ppv,"rU") as csv:
        firstline = csv.readline()
        raw_dict = defaultdict(list)
        refined_dict = defaultdict(list)

        for line in csv:

            line = line.strip().split(",")

            if line[2] != 'NA':

                ppv = float(line[2])

                if ppv == 0:
                    pass

                elif ppv > 0 and ppv <= 10:
                    if 10 not in raw_dict.keys():
                        to_add = [line[0],line[2]]
                        raw_dict[10] = list(to_add)
                    else:
                        raw_dict[10].append([line[0],line[2]])

                elif ppv > 10 and ppv <= 20:
                    if 20 not in raw_dict.keys():
                        to_add = [line[0],line[2]]
                        raw_dict[20] = list(to_add)
                    else:
                        raw_dict[20].append([line[0],line[2]])

                elif ppv > 20 and ppv <= 30:
                    if 30 not in raw_dict.keys():
                        to_add = [line[0],line[2]]
                        raw_dict[30] = list(to_add)
                    else:
                        raw_dict[30].append([line[0],line[2]])

                """elif ppv > 30 and ppv <= 40:
                    if 40 not in raw_dict.keys():
                        to_add = [line[0],line[2]]
                        raw_dict[40] = list(to_add)
                    else:
                        raw_dict[40].append([line[0],line[2]])"""


        max_key = max(raw_dict, key= lambda x: len(raw_dict[x]))
        max_values = raw_dict[max_key]

        q = 0

        while q < 2:

            if q == 0:
                trans = max_values[q]
                ppv = float(max_values[q+1])

                to_add = [trans,ppv]

                if ppv  > max_key-10 and ppv  <= max_key-8:
                    if max_key-8 not in refined_dict.keys():
                        raw_dict[max_key-8] = to_add
                    else:
                        raw_dict[max_key-8].append(to_add)

                elif ppv  > max_key-8 and ppv  <= max_key-6:
                    if max_key-6 not in refined_dict.keys():
                        refined_dict[max_key-6] = to_add
                    else:
                        refined_dict[max_key-6].append(to_add)

                elif ppv  > max_key-6 and ppv  <= max_key-4:
                    if max_key-4 not in refined_dict.keys():
                        refined_dict[max_key-4] = to_add
                    else:
                        refined_dict[max_key-4].append(to_add)

                elif ppv  > max_key-4 and ppv  <= max_key-2:
                    if max_key-2 not in refined_dict.keys():
                        refined_dict[max_key-2] = to_add
                    else:
                        refined_dict[max_key-2].append(to_add)

                elif ppv  > max_key-2 and ppv  <= max_key:
                    if max_key-2 not in refined_dict.keys():
                        refined_dict[max_key] = to_add
                    else:
                        refined_dict[max_key].append(to_add)

            elif q == 1:
                pass
            q+= 1

        for item in max_values[2:]:
            ppv = float(item[1])

            if ppv  > max_key-10 and ppv  <= max_key-8:
                if max_key-8 not in refined_dict.keys():
                    refined_dict[max_key-8] = item
                else:
                    refined_dict[max_key-8].append(item)

            elif ppv > max_key-8 and ppv <= max_key-6:
                if max_key-6 not in refined_dict.keys():
                    refined_dict[max_key-6] = [item[0]+item[1]]
                else:
                    refined_dict[max_key-6].append(item)

            elif ppv  > max_key-6 and ppv  <= max_key-4:
                if max_key-4 not in refined_dict.keys():
                    refined_dict[max_key-4] = item
                else:
                    refined_dict[max_key-4].append(item)

            elif ppv  > max_key-4 and ppv  <= max_key-2:
                if max_key-2 not in refined_dict.keys():
                    refined_dict[max_key-2] = item
                else:
                    refined_dict[max_key-2].append(item)

            elif ppv  > max_key-2 and ppv  <= max_key:
                if max_key not in refined_dict.keys():
                    refined_dict[max_key] = item
                else:
                    refined_dict[max_key].append(item)


        for key,value in refined_dict.iteritems():
            path = str(args.condition)+"_bin_"+str(key)

            try:
                os.mkdir(path,0755);
            except OSError:
                pass

            for item in value:
                if isinstance(item, float) == False:
                    if len(item) == 1:

                        contents = item.split(".")
                        transcript = contents[0]
                        command = 'cp ./PS/'+str(transcript)+"* "+path
                        subprocess.call(command,shell=True)

                    elif len(item) == 2:

                        contents = item[0].split(".")
                        transcript = contents[0]
                        command = 'cp ./PS/'+str(transcript)+"* "+path
                        subprocess.call(command,shell=True)

if __name__ == '__main__':
    main()
