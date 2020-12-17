import sys

with open(sys.argv[1],'r') as inp1:
    with open(sys.argv[2],'r') as inp2:
        with open(sys.argv[3],'w') as outp:
            input1 = inp1.readlines()
            input2 = inp2.readlines()

            shared = list(set(input1).intersection(input2))

            for item in shared:
                outp.write(item)
