import sys

transcript_list = []

with open(sys.argv[1], "rU") as handle:
    with open(sys.argv[2],"a+") as outp:

        firstline = handle.readline()
        outp.write(firstline)

        for record in handle:

            cont = record.split(",")

            sequence = cont[2]

            step = int(cont[4])

            if step > 0:
                step = step*10

            transcript = cont[3].split("|")

            new_contents = []

            for item in transcript:
                if ":" in str(item):
                    new_contents.append(item.split(":"))
                else:
                    new_contents.append(item)

            if 'CCAAT' in str(sequence) and 'ATG' in str(sequence) and int(new_contents[1][1]) - step < 31 and transcript not in transcript_list:
                #SD = sequence.index('CCAAT')

                #if sequence.index('ATG') > SD:
                outp.write(','.join(cont))
                transcript_list.append(transcript)
