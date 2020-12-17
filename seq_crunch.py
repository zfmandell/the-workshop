import sys

with open(sys.argv[1],"rU") as inp:
    firstline = inp.readline()
    potential = firstline.split(",")
    with open(sys.argv[2],"w") as outp:
        outp.write('transcript,FPUTR,CDS,TPUTR,')
        outp.write(",".join(potential[1:3]))
        outp.write(",sequence,")
        outp.write(",".join(potential[16:]))
        for line in inp:
            line = line.strip()
            line = line.split(",")

            sequence = line[3:16]
            seq = "".join(sequence)

            transcript_total = line[0].split("|")

            new_transcript = []

            for item in transcript_total:
                if ":" in str(item):
                    new_transcript.append(item.split(":"))
                else:
                    new_transcript.append(item)

            if str(new_transcript[1][0]) == 'FPUTR':
                FPUTR = new_transcript[1][1]
                CDS = new_transcript[2][1]
                if len(transcript_total) > 3:
                    TPUTR = new_transcript[3][1]
                else:
                    TPUTR = 'NA'

            elif str(new_transcript[1][0]) == 'CDS':
                FPUTR = 'NA'
                CDS = new_transcript[1][1]
                if len(transcript_total) > 2:
                    TPUTR = new_transcript[2][1]
                else:
                    TPUTR = 'NA'



            outp.write(str(new_transcript[0]))
            outp.write(","+str(FPUTR)+","+str(CDS)+","+str(TPUTR)+",")
            outp.write(",".join(line[1:3]))
            outp.write(","+str(seq)+",")
            outp.write(",".join(line[16:]))
            outp.write("\n")
