#!/usr/bin/bash


db="$1"
dbname="$2"
dbtype="$3"
file="$(readlink -f ${4})"
threads=$5
maxgb=$6

filedir="$(dirname ${file})"
filename="$(basename ${file})"
outdir="${filedir}/${dbname}"
mkdir -p "${outdir}"

outfile="${outdir}/${filename}.ANNOTATED.${dbtype}.fastq.gz"
echo "Annotating for : ${file} INTO ${outfile}"
bbduk.sh -"Xmx${maxgb}g" threads=$threads ref="${db}" in="${file}" int=f k=25 overwrite=t outm="${outfile}" maskmiddle=f rename=t
