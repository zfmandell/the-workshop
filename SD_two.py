import sys

with open(sys.argv[1], 'rU') as deaD:
    with open(sys.argv[2],'rU') as deaD_RMSD:
        with open(sys.argv[3],'rU') as deaD_deltaR:
            with open(sys.argv[4],"a+") as outp:

                firstline_deaD = deaD.readline()
                firstline_deaD = firstline_deaD.strip()
                outp.write(firstline_deaD+",nmrsd,deltaR")
                outp.write("\n")

                for element in deaD:
                    firstline_deaD_RMSD = deaD_RMSD.readline()
                    firstline_deaD_deltaR = deaD_deltaR.readline()

                    element = element.strip()

                    contents = element.split(",")

                    transcript = contents[0]

                    for item in deaD_RMSD:
                        item = item.strip()
                        contents_RMSD = item.split(",")

                        if str(contents_RMSD[0]) == str(transcript):

                            element = str(element)+","+str(contents_RMSD[1])


                    for item in deaD_deltaR:
                        item = item.strip()
                        contents_deltaR = item.split(",")

                        if str(contents_deltaR[0]) == str(transcript):

                            element = str(element)+","+str(contents_deltaR[1])

                    outp.write(element)
                    outp.write("\n")
                    deaD_RMSD.seek(0)
                    deaD_deltaR.seek(0)
