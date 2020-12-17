import sys

with open(sys.argv[1], 'rU') as DESeq:
    with open(sys.argv[2],'rU') as RMSD:
        with open(sys.argv[3],'rU') as deltaR:
            with open(sys.argv[4],"a+") as outp:

                firstline_DE = DESeq.readline()
                firstline_DE = firstline_DE.strip()
                outp.write("transcript,log2FoldChange,nmrsd,deltaR")
                outp.write("\n")

                for element in DESeq:

                    firstline_RMSD = RMSD.readline()
                    firstline_deltaR = deltaR.readline()

                    line = element.strip()
                    contents = line.split("\t")

                    transcript = contents[0]

                    transctipt_split = transcript.split("|")

                    transcript_final = []

                    for item in transctipt_split:
                        if ":" in str(item):
                            transcript_final.append(item.split(":"))
                        else:
                            transcript_final.append(item)

                    logchange = contents[3]

                    new_line = str(transcript_final[0])+","+str(logchange)

                    for item in RMSD:
                        item = item.strip()
                        contents_RMSD = item.split(",")

                        if str(contents_RMSD[0]) == str(transcript):

                            new_line = new_line+","+str(contents_RMSD[1])


                    for item in deltaR:
                        item = item.strip()
                        contents_deltaR = item.split(",")

                        if str(contents_deltaR[0]) == str(transcript):

                            new_line = new_line+","+str(contents_deltaR[1])

                    if len(new_line.split(",")) > 2:
                        outp.write(new_line)
                        outp.write("\n")

                    RMSD.seek(0)
                    deltaR.seek(0)
