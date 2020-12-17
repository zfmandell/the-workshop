#!/usr/bin/env python

from glob import glob
import os,sys
import subprocess
import shutil
import distutils.dir_util
from os import listdir
from os.path import isfile, join
import argparse


def is_int(value):
    try:
        int(value)
    except ValueError:
        return False
    except TypeError:
        return False

def find_bin(directory):
    q = 0
    while q < len(directory):
        sub_item = directory[q:q+2]
        if is_int(sub_item) != False:
            return sub_item
        q+= 1

def main():
    parser = argparse.ArgumentParser(description='merging ps of interest, MUST RUN IN FOLDER with ONLY directories from ppv_bin.py')
    parser.add_argument('control',default=None,help = '<.txt> biological control condition')
    parser.add_argument('experimental',default=None,help = '<.txt> biological experimental condition')

    args = parser.parse_args()

    directories = glob("./*/")

    directions = []

    for item in directories:
        directions.append(item.split("/")[1])

    total_files = []

    for item in directions:
        bin_num = find_bin(item)
        path = "./bin_"+str(bin_num)+"_"+str(args.control)+"_"+str(args.experimental)

        for direct in directories:
            if item in direct:
                onlyfiles = [f for f in listdir(direct) if isfile(join(direct, f))]

                if onlyfiles not in total_files:
                    total_files.append(onlyfiles)

        try:
            os.mkdir(path,0755);
        except OSError:
            pass

    done = []

    for item in directions:
        bin_num = find_bin(item)
        path = "./bin_"+str(bin_num)+"_"+str(args.control)+"_"+str(args.experimental)
        temp = []
        if bin_num not in done:
            for direct in directories:
                if bin_num in str(direct):
                    temp.append(direct)

            for files in total_files:
                for transcript in files:

                    tempfiles = [f for f in listdir(temp[0]) if isfile(join(temp[0], f))]

                    if transcript in tempfiles:

                        command_one = 'cat '+str(temp[0])+str(transcript) + " "+str(temp[1])+str(transcript)+">"+str(transcript)+"transcript_total.ps"
                        subprocess.call(command_one,shell=True)

                        command_two = "mv "+str(transcript)+"transcript_total.ps "+str(path)
                        subprocess.call(command_two,shell=True)

        done.append(bin_num)

if __name__ == '__main__':
    main()
