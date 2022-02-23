#!/usr/bin/env python3


import argparse
import numpy as np
from operator import itemgetter
from scipy.signal import savgol_filter

def consecutiveRanges(a, n):

	length = 1
	list = []
	
	# If the array is empty,
	# return the list
	if (n == 0):
		return list
	
	# Traverse the array
	# from first position
	for i in range (1, n + 1):
	
		# Check the difference
		# between the current
		# and the previous elements
		# If the difference doesn't
		# equal to 1 just increment
		# the length variable.
		if (i == n or a[i] -
			a[i - 1] != 1):
		
			# If the range contains
			# only one element.
			# add it into the list.
			if (length == 1):
				list.append(str(a[i - length]))
			else:
	
				# Build the range between the first
				# element of the range and the
				# current previous element as the
				# last range.
				temp = (str(a[i - length]) +
						"_" + str(a[i - 1]))
				list.append(temp)
		
			# After finding the
			# first range initialize
			# the length by 1 to
			# build the next range.
			length = 1
		
		else:
			length += 1
	return list

def read_coverage(cov_fyle):
	#create list of coverage values, position dependent
	return_list = []
	with open(cov_fyle,"r") as inp:
		for line in inp:
			return_list.append(float(line.strip().split('\t')[2]))
	return return_list

def transcript_check(first,second,cov,perc,len):
	if np.median(cov[first:second]) >= perc and second-first >= len:
		return True
	else:
		return False


def main():
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('fwd',type=str,help='')
	parser.add_argument('rev',type=str,help='')
	parser.add_argument('output',type=str,help='')
	parser.add_argument('-raw_threshold',type=int,default=5,help='raw cov threshold, default = 5')
	parser.add_argument('-threshold_percentile',type=int,default=25,help='percentile of peak coverage to be considered a transcript, default = 25')
	parser.add_argument('-threshold_length',type=int,default=500,help='minimum length to be considered a transcript, default = 500')
	parser.add_argument('-window_length',type=int,default=50,help='length of window when chunking, default = 50')
	parser.add_argument('-step_length',type=int,default=5,help='length of step when chunking, default = 5')
	
	args = parser.parse_args()

	Mutant_fwd = read_coverage(args.fwd)
	Mutant_fwd_smooth = list(savgol_filter(np.array(Mutant_fwd), 51, 3))
	Mutant_fwd_windows = [Mutant_fwd[i:i+args.window_length] for i in range(0, len(Mutant_fwd)-(args.window_length-1), args.step_length)]

	Mutant_perc_fwd = np.percentile(Mutant_fwd,args.threshold_percentile)

	fwd_keep = []
	for item in enumerate(Mutant_fwd_windows):
		if np.median(item[1]) >= args.raw_threshold:
			fwd_keep.append(item[0])

	fwd_transcripts_raw = consecutiveRanges(fwd_keep,len(fwd_keep))
	fwd_transcripts_refined = []
	for item in fwd_transcripts_raw:
		try:
			first = int(item.split('_')[0])*args.step_length
			second = int(item.split('_')[1])*args.step_length+args.window_length
			if transcript_check(first,second,Mutant_fwd,Mutant_perc_fwd,args.threshold_length) == True:
				fwd_transcripts_refined.append([first,second,'+'])
		except IndexError:
			pass

	Mutant_fwd,Mutant_fwd_windows,Mutant_fwd_smooth = [],[],[]

	Mutant_rev = read_coverage(args.rev)
	Mutant_rev_smooth = list(savgol_filter(np.array(Mutant_rev), 51, 3))
	Mutant_rev_windows = [Mutant_rev_smooth[i:i+args.window_length] for i in range(0, len(Mutant_rev_smooth)-(args.window_length-1), args.step_length)]

	Mutant_perc_rev = np.percentile(Mutant_rev,args.threshold_percentile)

	rev_keep = []
	for item in enumerate(Mutant_rev_windows):
		if np.median(item[1]) >= args.raw_threshold:
			rev_keep.append(item[0])

	rev_transcripts_raw = consecutiveRanges(rev_keep,len(rev_keep))
	rev_transcripts_refined = []
	for item in rev_transcripts_raw:
		try:
			first = int(item.split('_')[0])*args.step_length
			second = int(item.split('_')[1])*args.step_length+args.window_length
			if transcript_check(first,second,Mutant_rev,Mutant_perc_rev,args.threshold_length) == True:
				rev_transcripts_refined.append([first,second,'-'])
		except IndexError:
			pass

	Mutant_rev,Mutant_rev_windows,Mutant_rev_smooth = [],[],[]

	Mutant_final = fwd_transcripts_refined+rev_transcripts_refined
	Mutant_final_sorted = sorted(Mutant_final, key=itemgetter(0))

	with open(args.output,'w') as outp:
		outp.write('Strand,Start,End\n')
		
		for item in Mutant_final_sorted:
			Strand = str(item[2])
			if Strand == '+':
				Mutant_Start = str(item[0])
				Mutant_End = str(item[1])
			else:
				Mutant_Start = str(item[1])
				Mutant_End = str(item[0])

			outp.write(Strand+','+Mutant_Start+','+Mutant_End+'\n')

	
if __name__ == '__main__':
	main()







