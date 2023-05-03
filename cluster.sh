#!/bin/bash
#
dir=$2
clstr=1
cat $1 | awk -v dir=$2 -v str=1_contigs -v clstr=1 '{
        if (substr($0, 1, 1)==">") {filename=(substr($1,2) ".fasta")}
	print $0 >> dir"/"clstr"/"str"/"filename
	clstr=$clstr + 1
        close(filename)
}'

