#!/bin/bash
#
# Performs a local alignment on two sequences passed
# from the command line.
#
# Example:
#
#     local-align.sh THISLINE ISALIGNED -gapopen 1
#     local-align.sh THISLINE ISALIGNED -gapopen 10
#
# It also works if both entries are files.
#
#      local-align.sh file1.fa file2.fa  -gapopen 1
#
#
# Replace the tool 'water' with the tool named 'matcher'
# for a less rigorous but better performing algorithm.
# We use 'water' since it formats the output more nicely than 'matcher'.
#
ALIGNER=water

# Use this for higher performance.
# ALIGNER=matcher

# Run in strict mode.
set -ue

# Pass the remainder of the parameters at the with ${@:3}

if [ -f $1 ];
then
	$ALIGNER $1 $2 "${@:3}" -filter
else
	# Brad Pedersen: https://www.biostars.org/p/100676/#100689
	$ALIGNER <(echo -e ">a\n$1") <(echo -e ">b\n$2") "${@:3}" -filter
fi

