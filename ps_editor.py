#!/usr/bin/env python

import sys
from glob import glob
import argparse

def main():
    parser = argparse.ArgumentParser(description='binning and grabbing ps of interest, MUST RUN IN FOLDER WITH PS AND CT')
    parser.add_argument('condition',default=None,help = '<.txt> biological condition')

    args = parser.parse_args()

    ps = glob("./*.ps")

    cond = args.condition

    for item in ps:
        with open(item,"r+") as f:
            lines = f.readlines()

            for item in lines:
                insert = []
                contents = item.split(" ")

                if len(contents) == 9:
                    sub_item = contents[7].split(".")
                    if len(sub_item) > 1:
                        q = 0
                        index = lines.index(item)
                        while q < len(sub_item):
                            insert.append(sub_item[q])
                            if q == 0:
                                insert.append(cond)
                            q += 1
                        insert = ".".join(insert)

                        contents[7] = insert

                        lines[index] = contents



            f.seek(0)
            f.truncate()
            for item in lines:
                if not isinstance(item, basestring):
                    f.write(" ".join(item))
                else:
                    f.write(item)



if __name__ == '__main__':
    main()
