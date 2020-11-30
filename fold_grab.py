#!/usr/bin/env python


import sys
import glob
import subprocess

return_list = []
with open(sys.argv[1],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        return_list.append(line.strip())

for fyle in glob.glob(sys.argv[2]+'/*.ps'):
    if fyle.split("/")[4].split("_")[1] in return_list:
        subprocess.call(["cp",fyle,sys.argv[3]])
