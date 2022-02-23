#!/usr/bin/env python3

import sys
from collections import Counter

fasta = sys.argv[1]
output = sys.argv[2]

one,two,three,four = [],[],[],[]
with open(fasta,'r') as inp:
	for line in inp:
		if '>' in line:
			pass
		else:
			one.append(line.strip()[0])
			two.append(line.strip()[1])
			three.append(line.strip()[2])
			four.append(line.strip()[3])

one_c = Counter(one)
two_c = Counter(two)
three_c =Counter(three)
four_c = Counter(four)

one_f = [one_c['A']/len(one),one_c['C']/len(one),one_c['G']/len(one),one_c['U']/len(one)]
two_f = [two_c['A']/len(two),two_c['C']/len(two),two_c['G']/len(two),two_c['U']/len(two)]
three_f = [three_c['A']/len(three),three_c['C']/len(three),three_c['G']/len(three),three_c['U']/len(three)]
four_f = [four_c['A']/len(four),four_c['C']/len(four),four_c['G']/len(four),four_c['U']/len(four)]

with open(output,'w') as outp:
	outp.write('Pos,Nuc,Value\n')
	outp.write('1,A,'+str(one_f[0])+'\n')
	outp.write('1,C,'+str(one_f[1])+'\n')
	outp.write('1,G,'+str(one_f[2])+'\n')
	outp.write('1,U,'+str(one_f[3])+'\n')
	outp.write('2,A,'+str(two_f[0])+'\n')
	outp.write('2,C,'+str(two_f[1])+'\n')
	outp.write('2,G,'+str(two_f[2])+'\n')
	outp.write('2,U,'+str(two_f[3])+'\n')
	outp.write('3,A,'+str(three_f[0])+'\n')
	outp.write('3,C,'+str(three_f[1])+'\n')
	outp.write('3,G,'+str(three_f[2])+'\n')
	outp.write('3,U,'+str(three_f[3])+'\n')
	outp.write('4,A,'+str(four_f[0])+'\n')
	outp.write('4,C,'+str(four_f[1])+'\n')
	outp.write('4,G,'+str(four_f[2])+'\n')
	outp.write('4,U,'+str(four_f[3])+'\n')


