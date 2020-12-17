#!/usr/bin/env python

import glob
import os
import argparse
import subprocess



def main():
    parser = argparse.ArgumentParser(description='will dig into each folder, merge the two ps files corresponding to the same name')
    parser.add_argument('folds_1',type=str,help='name of first ps folder')
    parser.add_argument('folds_2',type=str,help='name of second ps folder')
    args = parser.parse_args()

    os.chdir(args.folds_1)
    names = glob.glob('*.ps')
    os.chdir("..")
    os.mkdir("final_ps",0755)
    for item in names:
        subprocess.call('cat '+str(args.folds_1)+"/"+str(item)+" "+str(args.folds_2)+"/"+str(item)+">"+"merged_"+str(item),shell=True)
    subprocess.call("mv *.ps final_ps",shell=True)

if __name__ == '__main__':
    main()
