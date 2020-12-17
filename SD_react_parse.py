#!/usr/bin/env python

import sys

with open(sys.argv[1],"r") as inp:
    with open("SD_"+str(sys.argv[1]),"a+") as outp:

        integer = 1

        prevline = ''

        for line in inp:
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

                if str(new_contents[1][0]) == 'FPUTR' and int(new_contents[1][1]) == 40:
                    base_contents = line.split("\t")

                    new_line = base_contents[:51]

                    header =  "|".join(contents)
                    outp.write(header)
                    outp.write('\n')
                    outp.write('\t'.join(new_line))
                    outp.write('\n')

                elif str(new_contents[1][0]) == 'FPUTR' and int(new_contents[1][1]) > 40:
                    base_contents = line.split("\t")

                    distance = int(new_contents[1][1])-40

                    new_line = base_contents[distance:distance+51]

                    header =  "|".join(contents)
                    outp.write(header)
                    outp.write('\n')
                    outp.write('\t'.join(new_line))
                    outp.write('\n')

                integer +=1
