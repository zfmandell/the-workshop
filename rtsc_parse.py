#!/usr/bin/env python

import os.path
import sys

react_input = raw_input("Please enter react file of interest: ")

if os.path.isfile(react_input) == False:
    print "File doesnt exist"
    sys.exit()

factor = False
while factor == False:
    react_output = raw_input("Please enter transcript feature: A) 5'UTR B) CDS C) 3'UTR D) 5'UTR+xCDS: ")
    nucleotides = ''
    if str(react_output.upper()) == 'A':
        choice = 'FPUTR'
        print "5'UTR selected"
        factor = True
    elif str(react_output.upper()) == 'B':
        choice = 'CDS'
        print "CDS selected"
        factor = True
    elif str(react_output.upper()) == 'C':
        choice = 'TPUTR'
        print "3'UTR selected"
        factor = True
    elif str(react_output.upper()) == 'D':
        choice = 'FPUTR_plus'
        nucleotides = raw_input("enter number of CDS nucleotides: ")
        print "5'UTR+"+nucleotides+"CDS selected"
        factor = True
    else:
        print "please choose A,B, or C"



with open(react_input,"r") as inp:
    with open(choice+nucleotides+"_"+react_input,"a+") as outp:

        integer = 1

        prevline = ''

        for line in inp:
            if line.strip() != '':
                line = line.strip()

                if integer % 2 != 0:
                    prevline = line
                    integer += 1

                elif integer % 2 == 0:

                    contents = prevline.split("|")

                    new_contents = []

                    for item in contents:
                        if ":" in str(item):
                            new_contents.append(item.split(":"))
                        else:
                            new_contents.append(item)

                    if choice == 'FPUTR' and str(new_contents[1][0]) == 'FPUTR':
                        base_contents = line.split("\t")

                        new_line = base_contents[0:int(new_contents[1][1])]

                        header = prevline
                        outp.write(header)
                        outp.write('\n')
                        outp.write('\t'.join(new_line))
                        outp.write('\n')

                    elif choice == 'FPUTR_plus' and str(new_contents[1][0]) == 'FPUTR':
                        base_contents = line.split("\t")

                        new_line = base_contents[0:int(new_contents[1][1])+int(nucleotides)]

                        header = prevline
                        outp.write(header)
                        outp.write('\n')
                        outp.write('\t'.join(new_line))
                        outp.write('\n')

                    elif choice == 'CDS' and str(new_contents[1][0]) == 'CDS':
                        base_contents = line.split("\t")

                        new_line = base_contents[0:int(new_contents[1][1])]

                        header = prevline
                        outp.write(header)
                        outp.write('\n')
                        outp.write('\t'.join(new_line))
                        outp.write('\n')

                    elif choice == 'CDS' and str(new_contents[2][0]) == 'CDS' and len(contents) < 4 :
                        base_contents=line.split("\t")

                        new_line = base_contents[int(new_contents[1][1]):]

                        header = prevline
                        outp.write(header)
                        outp.write('\n')
                        outp.write('\t'.join(new_line))
                        outp.write('\n')

                    elif choice == 'TPUTR' and len(contents) == 3:
                        if str(new_contents[2][0]) == 'TPUTR':
                            base_contents=line.split("\t")

                            new_line = base_contents[int(new_contents[1][1]):]

                            header = prevline
                            outp.write(header)
                            outp.write('\n')
                            outp.write('\t'.join(new_line))
                            outp.write('\n')

                    elif choice == 'TPUTR' and len(contents) > 3:
                        base_contents=line.split("\t")

                        number = (int(new_contents[1][1])) + (int(new_contents[2][1]))
                        new_line = base_contents[number:]

                        header = prevline
                        outp.write(header)
                        outp.write('\n')
                        outp.write('\t'.join(new_line))
                        outp.write('\n')

                    integer +=1
