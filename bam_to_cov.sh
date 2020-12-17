#!/bin/bash

set -e
set -u
set -o pipefail

find . -name "*.bam" > meta.txt

files=($(cut -f 1 meta.txt))

for bam_file in ${files[@]}
do
  #strip file extension, makes results
  results_file="$(basename $bam_file .bam).cov"

  #convert bam to cov
  bedtools genomecov -d  -ibam $bam_file -g ref.fai > ./cov/$results_file
done
