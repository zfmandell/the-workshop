#!/usr/bin/env python3

import sys

final = []
for line in sys.stdin:
	final.append(line.strip())

test = {}
for key,item in enumerate(final):
	if (key+3)%3 == 0:
		test[item.split(',')[0].strip('>')] = final[key+2].split(' ')[1].strip('(').strip(')')

print('coord,dG')
for key,value in test.items():
	print(str(key)+','+str(value))


	


		

  