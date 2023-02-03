#!/bin/bash
set -e
set -u

sample_id=${1}
reads=${2}

mkdir prefastqc_${sample_id}_logs
fastqc -o prefastqc_${sample_id}_logs --extract -f fastq -q ${reads}

# get sequence counts without using wc -l
grep "Total Sequences" prefastqc_${sample_id}_logs/*/*data.txt | cut -f 2 | awk -v prefix="${sample_id}\traw\t" '// { print prefix $0 }' > counts_${sample_id}_raw.tsv

