import sys

new_contents_1 = []
new_contents_2 = []

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


shared = list(set(new_contents_1).intersection(new_contents_2))

print shared
print len(shared)
