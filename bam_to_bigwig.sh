#!/bin/bash

set -e
set -u
set -o pipefail

find . -name "*.bam" > meta.txt

files=($(cut -f 1 meta.txt))

for bam_file in ${files[@]}
do
  #strip file extension, makes results
  results_file="$(basename $bam_file .bam).bw"

  #convert bam to bigwig
  bamCoverage -b $bam_file -o ../bigwig/$results_file -bs 3 -p 2
done
