import argparse
import operator

def index_build(index):
    'builds index dict'
    return_dict = {}
    with open(index,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            return_dict[float(line.strip().split(",")[0])] = [line.strip().split(",")[1],float(line.strip().split(",")[3])]
    return return_dict

def peak_owner(WT,KO):
    'builds list of shared peaks'
    shared = []
    for WT_key in WT.keys():
        for KO_key in KO.keys():
            if WT_key-3 <= KO_key <= WT_key+3:
                shared.append([WT_key,KO_key])
    WT_shared = [item[0] for item in shared]
    KO_shared = [item[1] for item in shared]

    WT_only = [x for x in WT.keys() if x not in WT_shared]
    KO_only = [x for x in KO.keys() if x not in KO_shared]

    return shared,WT_only,KO_only

def final_dict_build(shared,WT_only,KO_only,WT_index,KO_index):
    'build final dict to print'
    return_dict = {}
    for item in shared:
        return_dict[item[0]] = [WT_index[item[0]][0],WT_index[item[0]][1],KO_index[item[1]][1]]
    for item in WT_only:
        return_dict[item] = [WT_index[item][0],WT_index[item][1],'-']
    for item in KO_only:
        return_dict[item] = [KO_index[item][0],'-',KO_index[item][1]]

    return sorted(return_dict.iteritems(), key=lambda (k,v): operator.itemgetter(0)(v))

def main():
    parser = argparse.ArgumentParser(description='creates merged index file')
    parser.add_argument('WT',type=str,help='<.csv> WT delta index')
    parser.add_argument('KO',type=str,help='<.csv> KO delta index')
    parser.add_argument('outfile',type=str,default=None,help='name of output files')
    args = parser.parse_args()

    WT_index,KO_index = index_build(args.WT),index_build(args.KO)
    shared,WT_only,KO_only = peak_owner(WT_index,KO_index)
    final = final_dict_build(shared,WT_only,KO_only,WT_index,KO_index)

    with open(args.outfile,'w') as outp:
        outp.write('coordinate,context,WT_delta_KO_delta\n')
        for item in final:
            outp.write(str(item[0])+","+str(item[1][0])+","+str(item[1][1])+","+str(item[1][2])+"\n")

if __name__ == '__main__':
    main()
