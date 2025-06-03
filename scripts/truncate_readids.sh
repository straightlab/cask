#!/usr/bin/bash


file_dir="$1"
#/oak/stanford/groups/astraigh/charseq2.0/data/KAW_miseq_2024-01-24/data/KAW1/split_chimeras/SE_merge_pear/long
dr="$2"
#dna
infile="${file_dir}/${dr}.fastq.gz"
#dna.fastq.gz
outdir="${file_dir}/trunc_id"
outfile="${outdir}/${dr}.fastq"

echo "converting ${infile} to ${outfile}"

mkdir -p "${outdir}" && echo "Directory '${outdir}' created."

zcat "${infile}" | awk '{print $1}' > "${outfile}"

gzip "${outfile}"