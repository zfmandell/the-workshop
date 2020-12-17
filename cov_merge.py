import sys

def read_coverage(cov_fyle):
    #create list of coverage values, position dependent
    return_list = []
    with open(cov_fyle,"r") as inp:
        for line in inp:
            return_list.append(int(line.strip().split('\t')[2]))
    return return_list

one = read_coverage(sys.argv[1])
two = read_coverage(sys.argv[2])

three = [x + y for x, y in zip(one, two)]

i = 1
genome = 'NC_000964.3'

with open(sys.argv[3],'w') as outp:
    for item in three:
        outp.write(genome+'\t'+str(i)+'\t'+str(item)+"\n")
        i+=1
