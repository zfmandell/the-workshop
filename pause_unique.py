#!/usr/bin/env python3

# finds overlapping peaks of 3 barebones differential pause tables
# overlapping up and down are two distinct sets of outputs

import sys
import operator

dA = sys.argv[1]
dG = sys.argv[2]
dAdG = sys.argv[3]

def read_pauses(fyle):
	up_up,up_down,down_up,down_down = {},{},{},{}
	with open(fyle,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			coord = int(line.strip().split(',')[0])
			strand = line.strip().split(',')[1]
			Log2FC = line.strip().split(',')[10]
			seq = line.strip().split(',')[12]

			try:
				Log2FC = float(Log2FC)
				if Log2FC <= -2:
					up_down[coord] = [strand,abs(Log2FC),seq[:76]]
					down_down[coord] = [strand,abs(Log2FC),seq[76:]]
				elif Log2FC >= 2:
					up_up[coord] = [strand,abs(Log2FC),seq[:76]]
					down_up[coord] = [strand,abs(Log2FC),seq[76:]]
			except ValueError:
				pass
	return up_up,up_down,down_up,down_down


def pause_clean(pauses):
	temp,marked,return_dict = [],[],{}
	for key,value in pauses.items():
		temp.append([key,value[1]])
	temp_sorted = sorted(temp, key=operator.itemgetter(1))[::-1]
	for item in [x[0] for x in temp_sorted]:
		for subitem in [x[0] for x in temp_sorted]:
			if item == subitem:
				pass
			else:
				if item-20 <= subitem <= item+20:
					marked.append(subitem)
	temp_final = [x[0] for x in temp_sorted if x[0] not in marked]
	for item in temp_final:
		return_dict[item] = pauses[item]
	return return_dict

up_up,up_down,down_up,down_down = read_pauses(dA)

up_up_c = pause_clean(up_up)
up_down_c = pause_clean(up_down)
down_up_c = pause_clean(down_up)
down_down_c = pause_clean(down_down)

with open('up_up.fa','w') as outp:
	for key,value in up_up_c.items():
		outp.write('>'+str(key)+'_'+value[0]+'\n'+value[2]+'\n')

with open('up_down.fa','w') as outp:
	for key,value in up_down_c.items():
		outp.write('>'+str(key)+'_'+value[0]+'\n'+value[2]+'\n')

with open('down_up.fa','w') as outp:
	for key,value in down_up_c.items():
		outp.write('>'+str(key)+'_'+value[0]+'\n'+value[2]+'\n')

with open('down_down.fa','w') as outp:
	for key,value in down_down_c.items():
		outp.write('>'+str(key)+'_'+value[0]+'\n'+value[2]+'\n')




