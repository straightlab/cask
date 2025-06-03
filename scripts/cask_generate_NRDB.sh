#!/usr/bin/bash
db="$1" #unq or ci2
rootfolder="$2"
max_mem="$3" #max mem in gb (ex 350)
max_threads="$4" #max num threads



echo "Root folder is ${rootfolder}"


export LC_ALL="C"


cd "${rootfolder}/kdb"

### Generate NR unq/ci2 kmer db
echo "Making NRDB for ${db}"
if [ "$db" == "unq" ]; then
    kmc -k25 -m"${max_mem}" -sm -t"${max_threads}" -fm -ci1 -cx1 ../fa/repeats.fa ./k25/NRDB."${db}" tmp/
elif [ "$db" == "ci2" ]; then
    kmc -k25 -m"${max_mem}" -sm -t"${max_threads}" -fm -ci2 -cs100000000 ../fa/repeats.fa ./k25/NRDB."${db}" tmp/
else
    echo "Invalid db value: ${db}. Must be 'unq' or 'ci2'."
    exit 1
fi

### Interset NRDB with complement (BKG) to only keep kmers that are not present in the bkg
echo "Subtracting background from NRDB for ${db}"
kmc_tools simple ./k25/NRDB."${db}" ./k25/by_type/reptype_BKG kmers_subtract ./k25.minusBKG/NRDB."${db}"


### Dump NR db
echo "Dumping NRDB for ${db}" 
kmc_tools transform ./k25.minusBKG/NRDB."${db}" dump /dev/stdout -ci1 | awk 'BEGIN{OFS="\n"}{print ">"NR","$2, $1}' > ./k25.minusBKG/NRDB."${db}".dump.fa


### Add the UID from the NRDB to the "type-annotated" kmer db
echo "Making ADB for ${db}"
bbduk.sh -Xmx"${max_mem}"g threads="${max_threads}" ref=./k25.minusBKG/NRDB."${db}".dump.fa in=./k25.minusBKG/RDB.fa k=25 overwrite=t outm=./k25.minusBKG/ADB."${db}".fa.gz maskmiddle=f rename=t

### Make LUT
echo "Making LUT for ${db}"
zcat ./k25.minusBKG/ADB."${db}".fa.gz | awk -F $'\t' 'BEGIN{OFS=FS}(NR%2==1){split($1,x,","); split($2,y,","); print y[1], substr(x[1],2), x[3], substr(y[2],1,length(y[2])-2);}' | sort -k1,1n -k2,2 | awk -F $'\t' 'BEGIN{OFS=FS; ORS=""; getline; x=$1; print $1,$4, $2","$3;}{if ($1==x){print ","$2","$3;} else {x=$1; print "\n"$1, $4, $2","$3;}}END{print "\n"}' >./k25.minusBKG/LUT."${db}".txt

echo "Done with ${db}"