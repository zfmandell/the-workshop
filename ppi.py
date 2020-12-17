import sys

cspA_list_0 = []
cspA_list_1 = []
deaD_list_0 = []
deaD_list_1 = []
ID_map = []

def ID_switch(ID_mapped,transcript):
    ret = ''
    for item in ID_mapped:
        if str(transcript) in item:
            ret = item[6]
    return ret

with open("ecoli.txt","rU") as inp_3:
    for line in inp_3:
        line = line.strip()
        ID_map.append(line.split("\t"))

with open(sys.argv[1],"rU") as inp:
    first_line = inp.readline()
    for item in inp:
        line = item.strip()
        contents = line.split()
        cspA_list_0.append(contents[0])

with open(sys.argv[1],"rU") as inp:
    first_line = inp.readline()
    for item in inp:
        line = item.strip()
        contents = line.split()
        cspA_list_1.append(contents[1])

with open(sys.argv[2],"rU") as inp:
    first_line = inp.readline()
    for item in inp:
        line = item.strip()
        contents = line.split()
        deaD_list_0.append(contents[0])

with open(sys.argv[2],"rU") as inp:
    first_line = inp.readline()
    for item in inp:
        line = item.strip()
        contents = line.split()
        deaD_list_1.append(contents[1][9:])

cspA_final = []
deaD_final = []

cspA_final_switch = []
deaD_final_switch = []


for item in cspA_list_0:
    if item not in cspA_final:
        item = item.split(":")
        cspA_final_switch.append(ID_switch(ID_map,item[1]))
        cspA_final.append(item[1])

for item in cspA_list_1:
    if item not in cspA_final:
        item = item.split(":")
        cspA_final_switch.append(ID_switch(ID_map,item[1]))
        cspA_final.append(item[1])

for item in deaD_list_0:
    if item not in deaD_final:
        item = item.split(":")
        deaD_final_switch.append(ID_switch(ID_map,item[1]))
        deaD_final.append(item[1])

for item in deaD_list_1:
    if item not in deaD_final:
        item = item.split(":")
        deaD_final_switch.append(ID_switch(ID_map,item[1]))
        deaD_final.append(item[1])

shared = list(set(cspA_final).intersection(deaD_final))
shared_switch = list(set(cspA_final_switch).intersection(deaD_final_switch))

print ','.join(shared_switch)



print "len cspA "+str(len(cspA_final))
print "len deaD "+str(len(deaD_final))
print "len shared "+str(len(shared))


print shared_switch
