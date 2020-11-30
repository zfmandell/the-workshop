#!/bin/bash

set -e
set -u
set -o pipefail

find . -name "*.bam" > meta.txt

files=($(cut -f 1 meta.txt))

for bam_file in ${files[@]}
do
  #strip file extension, makes results
  results_file="$(basename $bam_file .bam).gff"

  #convert bam to cov
  stringtie $bam_file -G NC_000964.gff -o $results_file -p 2 -v
done
