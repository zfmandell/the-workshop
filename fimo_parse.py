import sys

with open(sys.argv[1], "rU") as handle:
    with open("fimo_output.csv","a+") as outp:

        outp.write("start,len,total")
        outp.write("\n")

        first_line = handle.readline()

        for fimo in handle:

            line = fimo.strip()

            line_cont = line.split("\t")

            contents = str(line_cont[2]).split("|")

            new_contents = []

            for item in contents:
                if ":" in str(item):
                    new_contents.append(item.split(":"))
                else:
                    new_contents.append(item)

            length = new_contents[1][1]

            if len(new_contents) > 3:

                total = int(new_contents[1][1])+int(new_contents[2][1])+int(new_contents[3][1])
                print total

            elif len(new_contents) > 2 and len(new_contents) < 4 :

                total = int(new_contents[1][1])+int(new_contents[2][1])
                print total

            else:

                total = int(new_contents[1][1])
                print total


            start = line_cont[3]

            outp.write(start)
            outp.write(",")
            outp.write(length)
            outp.write(",")
            outp.write(str(total))
            outp.write("\n")
