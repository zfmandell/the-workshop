import sys

models = [sys.argv[1],sys.argv[2]]

with open(sys.argv[1],"r") as one:
    firstline_one = one.readline()
    with open(sys.argv[2],'r') as two:
        firstline_two = two.readline()
        with open("merged_pearson_models.csv","w") as outp:

            outp.write(str(firstline_one)+"\n")

            one_dict = {}
            two_dict = {}

            for item_two in two:
                item_two = item_two.strip()
                item_two = item_two.split(",")
                transcript_two = item_two[0]

                two_dict[transcript_two] = item_two[1:]

            for item_one in one:
                item_one = item_one.strip()
                item_one = item_one.split(",")
                transcript_one = item_one[0]

                one_dict[transcript_one] = item_one[1:]


            for key, value in one_dict.iteritems():

                if str(key) in two_dict.keys():

                    outp.write(key)
                    outp.write(",")
                    for item in value:
                        outp.write(item)
                        outp.write(",")
                    outp.write("\n")

                    outp.write(key)
                    outp.write(",")
                    for item in two_dict[key]:
                        outp.write(item)
                        outp.write(",")
                    outp.write("\n")
