import argparse

def term_reader(csv):
	return_dict = {}
	with open(csv,'r') as inp:
        firstline = inp.readline()
        for line in inp:
        	return_dict[int(line.strip().split(',')[0])] = line.strip().split(',')[1]
    return return_dict

def rnet_reader(bw):
	return_plus,return_minus = {},{}
	with open(bw,'r') as inp:
		for line in inp:
			


def main():
    parser = argparse.ArgumentParser(description='determines 3-OH end idendity by combining Term-seq and RNET-seq information')
    parser.add_argument('Term-seq',type=str,help='list of Term-seq called 3-OH ends: <.csv>, no default value, must be inputted')
    parser.add_argument('RNET-seq',type=str,help='list of RNET-seq called 3-OH ends: <.bw>, no default value, must be inputted')
   
    args = parser.parse_args()

   

if __name__ == '__main__':
    main()
