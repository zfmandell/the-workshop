#!/usr/bin/env python3

import argparse

def term_reader(csv):
    return_dict = {}
    with open(csv,'r') as inp:
        firstline = inp.readline()
        for line in inp: 
            coord = int(line.strip().split(',')[0])
            strand = line.strip().split(',')[1]
            dist = int(line.strip().split(',')[2])
            return_dict[coord] = [strand,dist]
    
    return return_dict

def AU_calc(hairpin):
	i = 0
	A = hairpin[0:4]
	U = hairpin[-4:][::-1]
	for index,item in enumerate(A):
		
		if (item == 'A' or item == 'G') and (U[index] == 'U'):
			i += 1
		else:
			break
		
		"""
			

		if (item == 'A') and (U[index] == 'U'):
			if index == 0:
				i+=(1*7)
			elif index == 1:
				i+=(1*3.5)
			elif index == 2:
				i+=(1*0.5)
			elif index == 3:
				i+=(1*0.1875)
		elif (item == 'G') and (U[index] == 'U'):
			if index == 0:
				i+=(1*10)
			elif index == 1:
				i+=(1*5)
			elif index == 2:
				i+=(1*1)
			elif index == 3:
				i+=(1*0.25)
		"""
	
			

	return i

def table_reader(terms,table):
	return_dict = {}
	with open(table,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			coord = int(line.strip().split(',')[0])
			strand = line.strip().split(',')[1]
			hairpin = line.strip().split(',')[4]

			try:

				if terms[coord][0] == strand:
					AU = AU_calc(hairpin)
					return_dict[coord] = [strand,terms[coord][1],AU]
			except KeyError:
				pass 
			
	return return_dict
					
def main():
    parser = argparse.ArgumentParser(description='determines 3-OH end idendity by combining Term-seq and RNET-seq information')
    parser.add_argument('ends',type=str,help='list of Term-seq called 3-OH ends, adjusted using RNET-seq info: <.csv>, no default value, must be inputted')
    parser.add_argument('Total_Table',type=str,help='Table S2 from NusG Termination Paper: <.csv>, no default value, must be inputted')

    parser.add_argument('output',type=str,help='name of output file')
   
    args = parser.parse_args()

    terms = term_reader(args.ends)
    AU = table_reader(terms,args.Total_Table)

    with open(args.output,'w') as outp:
        outp.write('coord,strand,dist,AU\n')
        for key,value in AU.items():
            outp.write(str(key)+','+str(value[0])+','+str(value[1])+','+str(value[2])+'\n')
    
if __name__ == '__main__':
    main()