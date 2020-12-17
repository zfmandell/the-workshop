import sys

new_contents_1 = []
ID_map = []

def ID_switch(ID_mapped,transcript):
    ret = ''
    for item in ID_mapped:
        if str(transcript) in item:
            ret = item[2]
    return ret

with open(sys.argv[1],"rU") as inp_1:

    for item in inp_1:
        line = item.strip()
        contents = line.split("|")
        new_contents_1.append(contents[0])
        print new_contents_1

with open("ecoli.txt","rU") as inp_3:
    for line in inp_3:
        line = line.strip()
        ID_map.append(line.split("\t"))

with open("cspA_16_pos_switch.txt","a+") as outp_1:
    for element in new_contents_1:
        if ID_switch(ID_map,element) != False:

            prot_ID = ID_switch(ID_map,element)

            outp_1.write(prot_ID)
            outp_1.write("\n")
