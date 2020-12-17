import sys

DMS_genes = []
Gly_genes = []

DMS_1 = sys.argv[1]
Gly = sys.argv[2]
temp = sys.argv[3]

with open(DMS_1,"r") as DMS:
    with open(Gly,"r") as Gly:
        for line in DMS:
            contents = line.split("|")
            DMS_genes.append(contents[0])

        for line in Gly:
            contents = line.split("|")
            Gly_genes.append(contents[0])

sort_DMS = sorted(DMS_genes)
sort_Gly = sorted(Gly_genes)

total_list = set(sort_DMS).intersection(sort_Gly)

print len(total_list)

with open(DMS_1,"r") as DMS:
    with open(str(temp)+"_overlap.txt","a+") as outp:
        for line in DMS:
            contents = line.split("|")

            if contents[0] in total_list:
                outp.write(line)
