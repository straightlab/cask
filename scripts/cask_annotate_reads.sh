#!/usr/bin/bash
rootfolder="$1"
db="$2" #unq or ci2
input_fastq="$3"
outdir_suf="$4" #cen_kmatch/k25.minusBKG
max_mem="$5" #max mem in gb (ex 350)
max_threads="$6" #max num threads


file="$(readlink -f ${input_fastq})"
filename="$(basename ${file})"
ref="${rootfolder}/kdb/k25.minusBKG/NRDB.${db}.dump.fa"

outpath_pref="$(dirname ${input_fastq})"
full_outdir="${outpath_pref}/${outdir_suf}"
mkdir -p "${full_outdir}"
outfile="${full_outdir}/${filename}.ANNOTATED.${db}.fastq.gz"

echo "Starting annotation. Saving annotated reads to ${outfile}"

bbduk.sh -"Xmx${max_mem}g" threads="${max_threads}" ref="${ref}" in="${file}" int=f k=25 overwrite=t outm="${outfile}" maskmiddle=f rename=t