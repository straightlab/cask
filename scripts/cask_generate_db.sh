#!/usr/bin/bash
annotfile="$1"
genome_fa="$2"
genome_chrlength="$3"

rootfolder="$(pwd)"

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
kmc -k25 -m200 -t24 -fm -ci1 -cs100000000 ../fa/by_type/reptype_BKG.fa  ./k25/by_type/reptype_BKG tmp/kmc_tmpdir_BKG/

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

### Generate NR unq kmer and ci2 kmer db
echo "Making NRDB"
kmc -k25 -m200 -sm -t24 -fm -ci1 -cx1 ../fa/repeats.fa ./k25/NRDB.unq tmp/
kmc -k25 -m200 -sm -t24 -fm -ci2 -cs100000000 ../fa/repeats.fa ./k25/NRDB.ci2 tmp/

### Interset NRDB with complement (BKG) to only keep kmers that are not present in the bkg
echo "Subtracting background from NRDB"
kmc_tools simple ./k25/NRDB.unq ./k25/by_type/reptype_BKG kmers_subtract ./k25.minusBKG/NRDB.unq

kmc_tools simple ./k25/NRDB.ci2 ./k25/by_type/reptype_BKG kmers_subtract ./k25.minusBKG/NRDB.ci2

### Dump NR dbs
echo "Dumping NRDBS"
kmc_tools transform ./k25.minusBKG/NRDB.unq dump /dev/stdout -ci1 | awk 'BEGIN{OFS="\n"}{print ">"NR","$2, $1}' > ./k25.minusBKG/NRDB.unq.dump.fa

kmc_tools transform ./k25.minusBKG/NRDB.ci2 dump /dev/stdout -ci1 | awk 'BEGIN{OFS="\n"}{print ">"NR","$2, $1}' > ./k25.minusBKG/NRDB.ci2.dump.fa

### Add the UID from the NRDB to the "type-annotated" kmer db
echo "Making ADB"
bbduk.sh -Xmx200g threads=24 ref=./k25.minusBKG/NRDB.ci2.dump.fa in=./k25.minusBKG/RDB.fa k=25 overwrite=t outm=./k25.minusBKG/ADB.ci2.fa.gz maskmiddle=f rename=t

bbduk.sh -Xmx200g threads=24 ref=./k25.minusBKG/NRDB.unq.dump.fa in=./k25.minusBKG/RDB.fa k=25 overwrite=t outm=./k25.minusBKG/ADB.unq.fa.gz maskmiddle=f rename=t

### Make LUT
echo "Making LUT"
zcat ./k25.minusBKG/ADB.ci2.fa.gz | awk -F $'\t' 'BEGIN{OFS=FS}(NR%2==1){split($1,x,","); split($2,y,","); print y[1], substr(x[1],2), x[3], substr(y[2],1,length(y[2])-2);}' | sort -k1,1n -k2,2 | awk -F $'\t' 'BEGIN{OFS=FS; ORS=""; getline; x=$1; print $1,$4, $2","$3;}{if ($1==x){print ","$2","$3;} else {x=$1; print "\n"$1, $4, $2","$3;}}END{print "\n"}' >./k25.minusBKG/LUT.ci2.txt

zcat ./k25.minusBKG/ADB.unq.fa.gz | awk -F $'\t' 'BEGIN{OFS=FS}(NR%2==1){split($1,x,","); split($2,y,","); print y[1], substr(x[1],2), x[3], substr(y[2],1,length(y[2])-2);}' | sort -k1,1n -k2,2 | awk -F $'\t' 'BEGIN{OFS=FS; ORS=""; getline; x=$1; print $1,$4, $2","$3;}{if ($1==x){print ","$2","$3;} else {x=$1; print "\n"$1, $4, $2","$3;}}END{print "\n"}' >./k25.minusBKG/LUT.unq.txt

echo "Done."