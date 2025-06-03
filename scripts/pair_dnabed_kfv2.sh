#!/usr/bin/bash

#edited to include all reads, regardless of mapq & output gene_name and type

repdir="$1"
#/scratch/groups/astraigh/kelsey/centrochar/charseq_reads/K562/SUVDKO/TNG1/cen_kmatch/k25.minusBKG
dnabed="$2"
#/oak/stanford/groups/astraigh/charseq2.0/data/TNG_miseq_2024-04-01/data/TNG1/pairs/gencondeV29_hg38/all/dna.bed.gz


#kmersin fields: readid reptype_ag_ci2 class_aG_ci2 ci2kcount_inread reptype_ag_unq class_aG_unq unqkcount_inread
#kmersin="/scratch/groups/astraigh/kelsey/centrochar/charseq_reads/${rep}/cen_kmatch/k25.minusBKG/${rep}.dna.fastq.gz.ANNOTATED.ci2_and_unq.ag.txt"
kmersin="${repdir}/dna.fastq.gz.ANNOTATED.ci2_and_unq.ag.txt"


#outfields: readid ENSG reptype_ag_ci2 class_aG_ci2 ci2kcount_inread reptype_ag_unq class_aG_unq unqkcount_inread
foute="${repdir}/dnaside.bed.exons.CASK_mergev2.txt"
fouti="${repdir}/dnaside.bed.introns.CASK_mergev2.txt"



join -j1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3,2.4,2.5,2.6,2.7 <(zcat "${dnabed}" | awk -F $'\t' 'BEGIN{OFS=FS}(substr($12,4,1)=="T"){print $4,$7,$14,$15,$17}' | sort -k1,1) <(sort -k1,1 "${kmersin}") >"${foute}"



join -j1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3,2.4,2.5,2.6,2.7 <(zcat "${dnabed}" | awk -F $'\t' 'BEGIN{OFS=FS}(substr($12,4,1)=="G"){print $4,$7,$14,$15,$17}' | sort -k1,1) <(sort -k1,1 "${kmersin}") >"${fouti}"

