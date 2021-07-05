#!/usr/bin/bash


dbname="$1"
file="$(readlink -f ${2})"


filedir="$(dirname ${file})"
filename="$(basename ${file})"
outdir="${filedir}/${dbname}"

infile1="${outdir}/${filename}.ANNOTATED.ci2.ag.txt"
infile2="${outdir}/${filename}.ANNOTATED.unq.ag.txt"

outfile="${outdir}/${filename}.ANNOTATED.ci2_and_unq.ag.txt"
outfiletally="${outdir}/${filename}.ANNOTATED.ci2_and_unq.tally.txt"

export LC_ALL="C"
echo "Making joined db file ${outfile}"
join -e 0 -a 1 -a 2 -o 0,1.2,1.3,1.4,2.2,2.3,2.4 -t $'\t' <(sort -k1,1 "${infile1}") <(sort -k1,1 "${infile2}") >"${outfile}"

echo "Making tally table ${outfiletally}"
bedtools groupby -i <(cut -f2- "${outfile}" | sort -k1,1 -k2,2 -k4,4 -k5,5) -g 1,2,4,5 -c 1,3,6 -o count,sum,sum | sort -k5,5nr >"${outfiletally}"



