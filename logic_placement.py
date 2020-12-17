import argparse
import itertools

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('intrinsic',type=str)
    parser.add_argument('delta',type=str)
    args = parser.parse_args()

    delta_min = float(args.delta)

    n = 7
    lst = list(itertools.product([0, 1], repeat=n))
    dct = dict((item,0) for item in lst)

    with open(args.intrinsic,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            to_return = [0,0,0,0,0,0,0]

            if  float(line.strip().split(",")[5]) >= delta_min:
                to_return[0] = 1

            if  float(line.strip().split(",")[7]) >= delta_min:
                to_return[1] = 1

            if  float(line.strip().split(",")[9]) >= delta_min:
                to_return[2] = 1

            if  float(line.strip().split(",")[11]) >= delta_min:
                to_return[3] = 1

            if  float(line.strip().split(",")[13]) >= delta_min:
                to_return[4] = 1

            if  float(line.strip().split(",")[15]) >= delta_min:
                to_return[5] = 1

            if  float(line.strip().split(",")[17]) >= delta_min:
                to_return[6] = 1

            dct[tuple(to_return)] += 1
    with open('logic_final.csv','w') as outp:
        for key,value in sorted(dct.iteritems()):
            if value != 0:
                outp.write(', '.join(str(x) for x in key))
                outp.write(",")
                outp.write(str(value))
                outp.write("\n")

if __name__ == '__main__':
    main()
