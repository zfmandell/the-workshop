import glob

def longest(l):
    if(not isinstance(l, list)): return(0)
    return(max([len(l),] + [len(subl) for subl in l if isinstance(subl, list)]+[longest(subl) for subl in l]))

total = {}
for fyle in glob.glob('*.csv'):
    C,U,to_add = [],[],{}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            U.append(str(line.strip().split(",")[5].count('U')))
            C.append(str(line.strip().split(",")[5].count('C')+line.strip().split(",")[5].count('G')+line.strip().split(",")[5].count('A')))
        to_add['U'] = U
        to_add['C'] = C
    total[fyle.strip().split("_")[1]] = to_add


with open('C_pref.txt','w') as outp:
    outp.write('SI,dA,dG,dR\n')
    q = 0
    while q <= longest([total['SI']['C'],total['dA']['C'],total['dG']['C'],total['dR']['C']])-1:
        if q+1 <= len(total['SI']['C']):
            outp.write(total['SI']['C'][q])
            outp.write(",")
        else:
            outp.write('NA')
            outp.write(",")
        if q+1 <= len(total['dA']['C']):
            outp.write(total['dA']['C'][q])
            outp.write(",")
        else:
            outp.write('NA')
            outp.write(",")
        if q+1 <= len(total['dG']['C']):
            outp.write(total['dG']['C'][q])
            outp.write(",")
        else:
            outp.write('NA')
            outp.write(",")
        if q+1 <= len(total['dR']['C']):
            outp.write(total['dR']['C'][q])
            outp.write("\n")
        else:
            outp.write('NA')
            outp.write("\n")
        q+=1

with open('U_pref.txt','w') as outp:
    outp.write('SI,dA,dG,dR\n')
    q = 0
    while q <= longest([total['SI']['U'],total['dA']['U'],total['dG']['U'],total['dR']['U']])-1:
        if q+1 <= len(total['SI']['U']):
            outp.write(total['SI']['U'][q])
            outp.write(",")
        else:
            outp.write('NA')
            outp.write(",")
        if q+1 <= len(total['dA']['U']):
            outp.write(total['dA']['U'][q])
            outp.write(",")
        else:
            outp.write('NA')
            outp.write(",")
        if q+1 <= len(total['dG']['U']):
            outp.write(total['dG']['U'][q])
            outp.write(",")
        else:
            outp.write('NA')
            outp.write(",")
        if q+1 <= len(total['dR']['U']):
            outp.write(total['dR']['U'][q])
            outp.write("\n")
        else:
            outp.write('NA')
            outp.write("\n")
        q+=1
