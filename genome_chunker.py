#!/usr/bin/env python3

def genome_yield(fasta_name):
	seq = ''
	with open(fasta_name) as inp:
		for line in inp:
			if line[0] != '>':
				seq = seq+str(line.strip())
	return seq


def main():
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('sequence',type=str,help='')
	parser.add_argument('window',type=int,help='')
	parser.add_argument('step',type=int,help='')
	parser.add_argument('output',type=str,help='')

	args = parser.parse_args()

	genome = genome_yield(args.sequence)
	window = args.length
	step = args.step

	chunked = [genome[i:i+window] for i in xrange(0, len(genome)-(window-1), step)]

	
if __name__ == '__main__':
	main()
