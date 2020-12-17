import sys

deaD_dict = {}
cspA_dict = {}

shared_transcripts = {}
notshared_transcripts = {}

shared_final = []
notshared_final = []
ID_map = []

def uniq(lst):
    last = object()
    for item in lst:
        if item == last:
            continue
        yield item
        last = item
def ID_switch(ID_mapped,transcript):
    ret = ''
    for item in ID_mapped:
        if str(transcript) in item:
            ret = item[2]
    return ret

with open(sys.argv[1], "rU") as deaD:
    firstline = deaD.readline()

    for line in deaD:

        information = []

        line = line.strip()
        contents = line.split(",")
        change = contents[1]
        transcript = contents[3]
        step = contents[4]

        information.append(change)
        information.append(step)

        if transcript not in deaD_dict.keys():
            deaD_dict[transcript] = []
            deaD_dict[transcript].append(information)
        else:
            deaD_dict[transcript].append(information)

with open(sys.argv[2], "rU") as cspA:
    firstline = cspA.readline()

    for line in cspA:

        information = []

        line = line.strip()
        contents = line.split(",")
        change = contents[1]
        transcript = contents[3]
        step = contents[4]

        information.append(change)
        information.append(step)

        if transcript not in cspA_dict.keys():
            cspA_dict[transcript] = []
            cspA_dict[transcript].append(information)
        else:
            cspA_dict[transcript].append(information)


for transcript, information in deaD_dict.iteritems():

    if len(information) < 2:

        change_deaD, step_deaD = map(float, information[0])

        if transcript in cspA_dict.keys():

            if len(cspA_dict[transcript]) < 2:

                change_cspA, step_cspA = map(float, cspA_dict[transcript][0])

                if step_deaD == step_cspA:

                    if transcript not in shared_transcripts.keys():
                        shared_transcripts[transcript] = []
                        shared_transcripts[transcript].append(step_deaD)
                    else:
                        shared_transcripts[transcript].append(step_deaD)

                else:

                    if transcript not in notshared_transcripts.keys():
                        notshared_transcripts[transcript] = []
                        notshared_transcripts[transcript].append(step_deaD)
                    else:
                        notshared_transcripts[transcript].append(step_deaD)

            else:
                for pair in cspA_dict[transcript]:
                    change_cspA, step_cspA = map(float,pair)

                    if step_deaD == step_cspA:
                        if transcript not in shared_transcripts.keys():
                            shared_transcripts[transcript] = []
                            shared_transcripts[transcript].append(step_deaD)
                        else:
                            shared_transcripts[transcript].append(step_deaD)

                    else:

                        if transcript not in notshared_transcripts.keys():
                            notshared_transcripts[transcript] = []
                            notshared_transcripts[transcript].append(step_deaD)
                        else:
                            notshared_transcripts[transcript].append(step_deaD)


    else:

        for pair in information:
            change_deaD, step_deaD = map(float, pair)

            if transcript in cspA_dict.keys():

                if len(cspA_dict[transcript]) < 2:

                    change_cspA, step_cspA = map(float, cspA_dict[transcript][0])

                    if step_deaD == step_cspA:

                        if transcript not in shared_transcripts.keys():
                            shared_transcripts[transcript] = []
                            shared_transcripts[transcript].append(step_deaD)
                        else:
                            shared_transcripts[transcript].append(step_deaD)

                    else:

                        if transcript not in notshared_transcripts.keys():
                            notshared_transcripts[transcript] = []
                            notshared_transcripts[transcript].append(step_deaD)
                        else:
                            notshared_transcripts[transcript].append(step_deaD)

                else:
                    for pair in cspA_dict[transcript]:
                        change_cspA, step_cspA = map(float,pair)

                        if step_deaD == step_cspA:

                            if transcript not in shared_transcripts.keys():
                                shared_transcripts[transcript] = []
                                shared_transcripts[transcript].append(step_deaD)
                            else:
                                shared_transcripts[transcript].append(step_deaD)

                        else:

                            if transcript not in notshared_transcripts.keys():
                                notshared_transcripts[transcript] = []
                                notshared_transcripts[transcript].append(step_deaD)
                            else:
                                notshared_transcripts[transcript].append(step_deaD)

for transcript, information in shared_transcripts.iteritems():
    new_list = [str(int(i)) for i in information]
    shared_transcripts[transcript] = sorted(new_list)

for transcript, information in notshared_transcripts.iteritems():
    new_list = [str(int(i)) for i in information]
    notshared_transcripts[transcript] = sorted(new_list)

for item in shared_transcripts.keys():

    contents = item.split("|")

    new_contents = []

    for item in contents:
        if ":" in str(item):
            new_contents.append(item.split(":"))
        else:
            new_contents.append(item)

    shared_final.append(new_contents[0])

for item in notshared_transcripts.keys():

    contents = item.split("|")

    new_contents = []

    for item in contents:
        if ":" in str(item):
            new_contents.append(item.split(":"))
        else:
            new_contents.append(item)

    notshared_final.append(new_contents[0])


with open("ecoli.txt","rU") as inp_3:
    for line in inp_3:
        line = line.strip()
        ID_map.append(line.split("\t"))

with open("shared_37_final.txt","a+") as outp:
    outp.write("transcript,identical_steps")
    outp.write("\n")
    for transcript, information in shared_transcripts.iteritems():
        outp.write(transcript)
        outp.write("\t")
        outp.write(' '.join(information))
        outp.write(",")
        outp.write("\n")
with open("notshared_37_final.txt","a+") as outp:
    outp.write("transcript,identical_steps")
    outp.write("\n")
    for transcript, information in notshared_transcripts.iteritems():
        outp.write(transcript)
        outp.write("\t")
        outp.write(' '.join(information))
        outp.write(",")
        outp.write("\n")

with open("shared_37.txt","a+") as outp_1:

    for element in shared_final:

        if ID_switch(ID_map,element) != False:

            prot_ID = ID_switch(ID_map,element)

            outp_1.write(prot_ID)
            outp_1.write("\n")

with open("notshared_37.txt","a+") as outp_2:

    for element in notshared_final:
        if ID_switch(ID_map,element) != False:

            prot_ID = ID_switch(ID_map,element)

            outp_2.write(prot_ID)
            outp_2.write("\n")

print
