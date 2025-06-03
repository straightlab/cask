#!/usr/bin/bash
annotfile="$1"
genome_fa="$2"
genome_chrlength="$3"
rootfolder="$4"
max_mem="$5" #max mem in gb ie 350


echo "Root folder is ${rootfolder}"

echo "Extracting repeat classes"

export LC_ALL="C"
sort -k1,1 -k2,2n "$annotfile" >annotations.bed

awk -F $'\t' 'BEGIN{OFS=FS}{print $1, $2, $3, $4, $5, $6}' annotations.bed >annotations.forfa.bed 

awk -F $'\t' 'BEGIN{OFS=FS}{print $4, $7}' annotations.bed | sort -k1,1 -k2,2 | uniq >repeats.withClass.txt
cut -f1,1 repeats.withClass.txt | sort -k1,1 | uniq >repeats.txt


### Make fa
echo "Making fastas"
mkdir -p fa

bedtools getfasta -nameOnly -fo fa/repeats.fa -fi "${genome_fa}" -bed annotations.forfa.bed


### Split fa
echo "Splitting fastas"
cd "${rootfolder}/fa"
awk_body='BEGIN{OFS=","}(NR==FNR){x[NR]=$1; y[$1]=NR; k[$1]=0; next;}(/^>/){z=substr($1,2); if (z in y){fname="by_type/reptype_"z".fa"; print ">"z, k[z], FNR > fname; getline; print $0 > fname; k[z]+=1;}}'
echo $awk_body >awkbody.awk

mkdir by_type
awk -f awkbody.awk ../repeats.txt repeats.fa

### Generate kmer db for each repeat type
echo "Generating kmer dbs"


cd "${rootfolder}"
mkdir -p ./kdb/k25/by_type

cd "${rootfolder}/kdb"

parallel -t -j1 'mkdir -p tmp/kmc_tmpdir_{}/' :::: ../repeats.txt

parallel -t -j12 'kmc -k25 -m12 -t2 -fm -ci1 -cs100000000 ../fa/by_type/reptype_{}.fa  ./k25/by_type/reptype_{} tmp/kmc_tmpdir_{}/' :::: ../repeats.txt


### Generate kmer db for the complement of the annotation file
echo "Generating background db"
bedtools complement -i ../annotations.forfa.bed -g "${genome_chrlength}" | sort -k1,1 -k2,2n | awk -F $'\t' 'BEGIN{OFS=FS}{print $1, $2, $3, "BKG"}' >../annotations.BKG.bed

bedtools getfasta -nameOnly -fo ../fa/by_type/reptype_BKG.fa -fi "${genome_fa}" -bed ../annotations.BKG.bed


mkdir -p tmp/kmc_tmpdir_BKG
kmc -k25 -m"${max_mem}" -t24 -fm -ci1 -cs100000000 ../fa/by_type/reptype_BKG.fa  ./k25/by_type/reptype_BKG tmp/kmc_tmpdir_BKG/

### Intersec the by_type kmer db with complement(BKG) to only keep kmers that are not present in the bkg
echo "Subtracting background"
mkdir -p ./k25.minusBKG/by_type/
parallel -j24 -t 'kmc_tools simple ./k25/by_type/reptype_{} ./k25/by_type/reptype_BKG kmers_subtract ./k25.minusBKG/by_type/reptype_{}' :::: ../repeats.txt

### Dump the databases into fa files
echo "Dumping type dbs to fa"
mkdir -p ./k25.minusBKG/by_type/dump/

parallel -j24 -t 'kmc_tools transform ./k25.minusBKG/by_type/reptype_{} dump ./k25.minusBKG/by_type/dump/reptype_{}.txt -ci1' :::: ../repeats.txt

awkbody='BEGIN{OFS="\n"}{print ">"reptype","NR","$2, $1}'
echo $awkbody >awkbody.awk

parallel -j24 -t 'awk -v reptype="{}" -f awkbody.awk ./k25.minusBKG/by_type/dump/reptype_{}.txt >./k25.minusBKG/by_type/dump/reptype_{}.fa' :::: ../repeats.txt


### Combine all the fasta dumped databases into a single fa file
echo "Making RDB"
touch ./k25.minusBKG/RDB.fa
rm ./k25.minusBKG/RDB.fa
touch ./k25.minusBKG/RDB.fa
parallel -j1 -t 'cat ./k25.minusBKG/by_type/dump/reptype_{}.fa >>./k25.minusBKG/RDB.fa' :::: ../repeats.txt 


echo "Done with RDB"