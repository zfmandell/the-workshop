#!/usr/bin/env python

from __future__ import division
from glob import glob
import argparse

def main():
    parser = argparse.ArgumentParser(description='generate csv of stats regarding intra/inter feature base pairing for every CT in a file')
    parser.add_argument('name',default=None,help = '<.txt> intended name of csv file')

    args = parser.parse_args()

    files = glob("./*.ct")

    with open(args.name,"a+") as outp:
        outp.write("transcript,FPUTR_bp,CDS_bp,TPUTR_bp\n")

    for f in files:
        with open(f,"r") as inp:
            firstline = inp.readline()
            firstline = firstline.strip()
            contents = firstline.split()
            trans_cont = contents[1].split(".")
            new_trans_cont = []

            for item in trans_cont:
                if ":" in str(item):
                    new_trans_cont.append(item.split(":"))
                else:
                    new_trans_cont.append(item)

            if str(new_trans_cont[1][0]) == 'FPUTR':
                FPUTR = int(new_trans_cont[1][1])
            else:
                FPUTR = 0
            if str(new_trans_cont[1][0]) == 'CDS':
                CDS = int(new_trans_cont[1][1])
            elif str(new_trans_cont[2][0]) == 'CDS':
                CDS = int(new_trans_cont[2][1])
            if len(trans_cont) == 3:
                if str(new_trans_cont[2][0]) == 'TPUTR':
                    TPUTR = int(new_trans_cont[2][1])
                else:
                    TPUTR = 0
            elif len(trans_cont) == 4:
                TPUTR = int(new_trans_cont[3][1])

            FPUTR_total = 0
            CDS_total = 0
            TPUTR_total = 0
            line_count = 1

            for line in inp:
                line = line.strip()
                info = line.split()

                if FPUTR != 0 and line_count <= FPUTR:
                    if int(info[4]) > 0:
                        FPUTR_total += 1
                elif FPUTR != 0 and line_count > FPUTR and len(trans_cont) == 3:
                    if int(info[4]) > 0:
                        CDS_total += 1
                elif FPUTR != 0 and line_count > FPUTR and line_count <= FPUTR+CDS and len(trans_cont) == 4:
                    if int(info[4]) > 0:
                        CDS_total += 1
                elif FPUTR != 0 and line_count > FPUTR and line_count > FPUTR+CDS and len(trans_cont) == 4:
                    if int(info[4]) > 0:
                        TPUTR_total += 1
                elif FPUTR == 0:
                    if line_count <= CDS and TPUTR == 0:
                        if int(info[4]) > 0:
                            CDS_total += 1
                    elif line_count <= CDS and TPUTR != 0:
                        if int(info[4]) > 0:
                            CDS_total += 1
                    elif line_count > CDS and TPUTR != 0:
                        if int(info[4]) > 0:
                            TPUTR_total = 0

                line_count += 1

        with open(args.name,"a+") as outp:
            if FPUTR == 0 and TPUTR == 0 :
                FPUTR_write = 'NA'
                CDS_write = CDS_total
                TPUTR_write = 'NA'
            elif FPUTR == 0 and TPUTR != 0:
                FPUTR_write = 'NA'
                CDS_write = CDS_total
                TPUTR_write = TPUTR_total
            elif FPUTR != 0 and TPUTR == 0:
                FPUTR_write = FPUTR_total
                CDS_write = CDS_total
                TPUTR_write = 'NA'
            else:
                FPUTR_write = FPUTR_total
                CDS_write = CDS_total
                TPUTR_write = TPUTR_total

            feature_list = [FPUTR_write,CDS_write,TPUTR_write]

            outp.write(str(".".join(trans_cont)))
            outp.write(",")

            for item in feature_list:
                outp.write(str(item))
                outp.write(",")
            outp.write("\n")


if __name__ == '__main__':
    main()
