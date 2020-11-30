#!/bin/bash
#
# Performs a local alignment on two sequences passed
# from the command line.
#
# Example:
#
#     global-align.sh THISLINE ISALIGNED -gapopen 1
#     global-align.sh THISLINE ISALIGNED -gapopen 10
#
# It also works if both entries are files.
#
#      global-align.sh file1.fa file2.fa  -gapopen 1
#
#
# Replace the tool 'needle' with the tool named 'stretcher'
# for a less rigorous but better performing algorithm.
# We use needle since it formats the output more nicely than strecher.
#

ALIGNER=needle

# Uncomment this for higher performance.
# ALIGNER=strecher

# Run in strict mode.
set -ue

# Pass the remainder of the parameters at the with ${@:3}

if [ -f $1 ];
then
	needle $1 $2 "${@:3}" -filter
else
	# Brad Pedersen: https://www.biostars.org/p/100676/#100689
	needle <(echo -e ">a\n$1") <(echo -e ">b\n$2") "${@:3}" -filter
fi

