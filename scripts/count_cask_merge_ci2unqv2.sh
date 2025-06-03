#!/usr/bin/bash

export LC_COLLATE="C"

infile="$1"

outfile_prefix="${infile%.txt}"
echo "Input: $infile"

reptype_outf="${outfile_prefix}.reptypecount_ci2_and_unqv2.txt"
class_outf="${outfile_prefix}.classcount_ci2_and_unqv2.txt"




#reptype
cut -d $'\t' -f 2,3,4,5,6,9 "$infile" | sort | uniq -c | sed 's/^\ *//' | sed 's/\ /\t/' | sort -k2,2 -k1,1nr > "$reptype_outf"
echo "Output 1/2: $reptype_outf"
#class
cut -d $'\t' -f 2,3,4,5,7,10 "$infile" | sort | uniq -c | sed 's/^\ *//' | sed 's/\ /\t/' | sort -k2,2 -k1,1nr > "$class_outf"

echo "Output 2/2: $class_outf"
