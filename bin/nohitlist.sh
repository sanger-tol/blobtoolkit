#!/bin/bash

# input
fasta=$1
blast=$2
prefix=$3
E=$4

# find ids of sequences with no hits in the blastx search
grep '>' $fasta | \
    grep -v -w -f <(awk -v evalue="$E" '{{if($14<{evalue}){{print $1}}}}' $blast | sort | uniq) | \
    cut -f1 | sed 's/>//' > $prefix.nohit.txt



