#!/usr/bin/bash


repdir="$1"
#/scratch/groups/astraigh/kelsey/centrochar/charseq_reads/K562/SUVDKO/TNG1/cen_kmatch/k25.minusBKG
dnabed="$2"
#/oak/stanford/groups/astraigh/charseq2.0/data/TNG_miseq_2024-04-01/data/TNG1/pairs/gencondeV29_hg38/all/dna.bed.gz


foute="${repdir}/dnaside.bed.exons.countallRNAs.txt"
fouti="${repdir}/dnaside.bed.introns.countallRNAs.txt"

zcat "$dnabed" > "${repdir}/allpairs.dnaside.bed.tmp"

wait


awk -F $'\t' 'BEGIN{OFS=FS}(substr($12,4,1)=="T"){print $17}' "${repdir}/allpairs.dnaside.bed.tmp" | sort -k1,1 |  uniq -c | sed 's/^\ *//' | sed 's/\ /\t/' | sort -k2,2 -k1,1nr>"${foute}" &


awk -F $'\t' 'BEGIN{OFS=FS}(substr($12,4,1)=="G"){print $17}' "${repdir}/allpairs.dnaside.bed.tmp" | sort -k1,1 |  uniq -c | sed 's/^\ *//' | sed 's/\ /\t/' | sort -k2,2 -k1,1nr>"${fouti}" &

wait

rm "${repdir}/allpairs.dnaside.bed.tmp"