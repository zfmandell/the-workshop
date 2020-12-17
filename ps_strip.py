
import glob
import sys

ends,seq_dict = [],{}

for item in glob.glob("./not_sure/*.ps"):
    with open(item,"r") as f:
        for item in f:
            contents = item.split(" ")
            if len(contents) == 9:
                sub_item = contents[7].split(".")
                if len(sub_item) > 1:
                    ends.append(int(sub_item[0].split("_")[0].split(":")[1]))

for item in glob.glob("./CT/*.ct"):
    with open(item,"r") as f:
        firstline = f.readline()
        temp = []
        for line in f:
            temp.append(line.split()[1])
        print len(temp)
        seq_dict[int(firstline.split()[1].split("_")[0].split(":")[1])] = ''.join(str(e) for e in temp[:-1])

with open(sys.argv[1],'w') as outp:
    outp.write("POT,Seq\n")
    for key,value in seq_dict.iteritems():
        if key in ends:
            outp.write(str(key)+','+str(value)+"\n")
