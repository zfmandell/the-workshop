#!/usr/bin/env python

import argparse


def main():
    parser = argparse.ArgumentParser(description='Generate merged overlap file')
    parser.add_argument("overlap_1",type=str,help="overlap file 1")
    parser.add_argument("overlap_2",type=str,help="overlap file 2")
    parser.add_argument('name',type=str, help='name of new overlap file')

    args = parser.parse_args()

    with open(args.overlap_1,"rU") as over_1:
        over_1 = [line.strip() for line in over_1]

    with open(args.overlap_2,"rU") as over_2:
        over_2 = [line.strip() for line in over_2]

    shared = list(set(over_1).intersection(over_2))

    with open(args.name,"w") as outp:
        for item in shared:
            outp.write(str(item)+"\n")

if __name__ == '__main__':
    main()
