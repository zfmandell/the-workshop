import sys

new_contents_1 = []
new_contents_2 = []
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

with open(sys.argv[2],"rU") as inp_2:

    for item in inp_2:
        line = item.strip()
        contents = line.split("|")
        new_contents_2.append(contents[0])

with open("ecoli.txt","rU") as inp_3:
    for line in inp_3:
        line = line.strip()
        ID_map.append(line.split("\t"))

with open("deaD_16.txt","a+") as outp_1:

    for element in new_contents_1:

        if element not in new_contents_2:

            if ID_switch(ID_map,element) != False:

                prot_ID = ID_switch(ID_map,element)

                outp_1.write(prot_ID)
                outp_1.write("\n")

with open("cspA_16.txt","a+") as outp_2:
    for element in new_contents_2:


        if element not in new_contents_1:

            if ID_switch(ID_map,element) != '':

                prot_ID = ID_switch(ID_map,element)

                outp_2.write(prot_ID)
                outp_2.write("\n")


with open("shared_16.txt","a+") as outp_3:
    shared = list(set(new_contents_1).intersection(new_contents_2))

    for element in shared:

        if ID_switch(ID_map,element) != False:

            prot_ID = ID_switch(ID_map,element)

            outp_3.write(prot_ID)
            outp_3.write("\n")
