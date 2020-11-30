import argparse
import glob

def main():
    parser = argparse.ArgumentParser(description='compile MFE from different 5prime lengths')
    parser.add_argument('index',type=str,default=None,help='<.txt> index')
    parser.add_argument('outname',type=str,default=None,help='output .csv name')
    args = parser.parse_args()

    ends = {}
    with open(args.index,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            ends[float(line.strip().split(",")[0])] = [line.strip().split(",")[1],float(line.strip().split(",")[3])]

    all_len = []
    for fyle in glob.glob('*MFE.csv*')[1:]:
        size = fyle.split("_")[0]
        print size
        size = {}
        with open(fyle,'r') as inp:
            firstline = inp.readline()
            for line in inp:
                size[float(line.strip().split(",")[0].split(":")[1])] = float(line.strip().split(",")[1])
            all_len.append(size)
    with open(glob.glob('*MFE.csv*')[0],'r') as inp:
        firstline = inp.readline()
        size = glob.glob('*MFE.csv*')[0].split("_")[0]
        print size
        size = {}
        for line in inp:
            size[float(line.strip().split(",")[0].split(":")[1])] = float(line.strip().split(",")[1])
        all_len.append(size)

    with open(args.outname,'w') as outp:
        outp.write('end_coordinate,feature,abundance,deltaG_20,deltaG_30,deltaG_40,deltaG_50,deltaG_60,deltaG_70,deltaG_80,deltaG_90,deltaG_100\n')
        for key,value in ends.iteritems():
            if value[1] > 100:
                outp.write(str(key)+","+str(value[0])+","+str(value[1]))
                for item in all_len:
                    outp.write(","+str(item[key]))
                outp.write("\n")


if __name__ == '__main__':
    main()
