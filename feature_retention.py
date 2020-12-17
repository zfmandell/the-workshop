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
        outp.write("transcript,intra_FPUTR,inter_FPUTR,intra_CDS,inter_CDS,intra_TPUTR,inter_TPUT,total_bp_FPUTR,total_bp_CDS,total_bp_FPUTR\n")

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

            FPUTR_intra = 0
            CDS_intra = 0
            TPUTR_intra = 0
            FPUTR_inter = 0
            CDS_inter= 0
            TPUTR_inter = 0
            line_count = 1

            for line in inp:
                line = line.strip()
                info = line.split()

                if FPUTR != 0 and line_count <= FPUTR:
                    if int(info[4]) <= FPUTR and int(info[4]) > 0:
                        FPUTR_intra += 1
                    elif int(info[4]) > FPUTR and int(info[4]) > 0:
                        FPUTR_inter += 1
                elif FPUTR != 0 and line_count > FPUTR and len(trans_cont) == 3:
                    if int(info[4]) > FPUTR:
                        CDS_intra += 1
                    else:
                        CDS_inter += 1
                elif FPUTR != 0 and line_count > FPUTR and line_count <= FPUTR+CDS and len(trans_cont) == 4:
                    if int(info[4]) > FPUTR and int(info[4]) <= FPUTR+CDS:
                        CDS_intra += 1
                    elif int(info[4]) <= FPUTR or int(info[4]) > FPUTR+CDS:
                        CDS_inter += 1
                elif FPUTR != 0 and line_count > FPUTR and line_count > FPUTR+CDS and len(trans_cont) == 4:
                    if int(info[4]) > FPUTR+CDS:
                        TPUTR_intra += 1
                    else:
                        TPUTR_inter += 1
                elif FPUTR == 0:
                    if line_count <= CDS and TPUTR == 0:
                        CDS_intra += 1
                    elif line_count <= CDS and TPUTR != 0:
                        if int(info[4]) <= CDs:
                            CDS_intra += 1
                        else:
                            CDS_inter += 1
                    elif line_count > CDS and TPUTR != 0:
                        if int(info[4]) > CDs:
                            TPUTR_intra = 0
                        else:
                            TPUTR_inter +=1

                line_count += 1

        with open(args.name,"a+") as outp:
            if FPUTR == 0 and TPUTR == 0 :
                FPUTR_intra = 'NA'
                FPUTR_inter = 'NA'
                TPUTR_intra = 'NA'
                TPUTR_inter = 'NA'
                total_bp_FPUTR = 'NA'
                total_bp_CDS = CDS_inter+CDS_intra
                total_bp_TPUTR = 'NA'
            elif FPUTR == 0 and TPUTR != 0:
                FPUTR_intra = 'NA'
                FPUTR_inter = 'NA'
                total_bp_FPUTR = 'NA'
                total_bp_CDS = CDS_inter+CDS_intra
                total_bp_TPUTR = TPUTR_inter+TPUTR_intra
            elif FPUTR != 0 and TPUTR == 0:
                TPUTR_intra = 'NA'
                TPUTR_inter = 'NA'
                total_bp_FPUTR = FPUTR_inter+FPUTR_intra
                total_bp_CDS = CDS_inter+CDS_intra
                total_bp_TPUTR = 'NA'
            else:
                total_bp_FPUTR = FPUTR_inter+FPUTR_intra
                total_bp_CDS = CDS_inter+CDS_intra
                total_bp_TPUTR = TPUTR_inter+TPUTR_intra

            feature_list = [FPUTR_intra,FPUTR_inter,CDS_intra,CDS_inter,TPUTR_intra,TPUTR_inter,total_bp_FPUTR,total_bp_CDS,total_bp_TPUTR]

            outp.write(str(".".join(trans_cont)))
            outp.write(",")

            for item in feature_list:
                outp.write(str(item))
                outp.write(",")
            outp.write("\n")





if __name__ == '__main__':
    main()
