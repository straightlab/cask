# CASK
Classification of Ambivalent Sequences using K-mers.

CASK is a tool to classify sequencing reads from repetitive genomic regions based on their k-mer composition. CASK works by first generating an annotated k-mer database, linking all the k-mers in the genome with an ambivalence group, which is the set of all repeat types in which the k-mer occurs. CASK then uses the annotated k-mer database to parse sequencing reads and annotate them with an ambivalence group consistent with all the kmers in the read. 

CASK relies on KMC for generating repeat-specific k-mer databases and on BBduk for extracting relevant kmers from the reads. Generation of the annotated k-mer databases and annotation of the reads is currently done through bash scripts and a small python module.

The inputs for CASK are :
- a fasta file with the genome sequence
- a bed file with the genomic coordinates of the repeats and their type (obtained typically with RepeatMasker)
- a fastq file with sequencing reads to be classified

## Dependencies
Set up and activate conda environment:  
`conda env create -f cask_env.yml`
`conda activate cask`

OR the following tools need to be installed separately.
- bedtools (https://bedtools.readthedocs.io/en/latest/)
- KMC (https://github.com/refresh-bio/KMC)
- bbduk.sh from the BBMap suite (https://sourceforge.net/projects/bbmap/) 


# CASK run instructions


## Step 0: Set up working directories & copy over files

1. Output directory for k-mer databases (independent of reads)
    
    Create a directory to keep the k-mer databases. You can use these same k-mer databases to annotate multiple sequence datasets, so pick a location independent of the reads. For example: `/scratch/groups/astraigh/kelsey/centrochar/cen_only`

2. Working directory for read annotation and downstream analysis
    
    Create a directory where you want to store your annotated reads files as well as all the downstream analysis files. This can either be in the same directory as the raw reads or you can create a symlink to your reads in a different directory. Downstream steps will create more directories within this so just create the rootdir. For example:
    `/scratch/groups/astraigh/kelsey/centrochar/charseq_reads/decon/WT_K562/WTrep1`
3. Copy or create symlink to raw reads
    
    `ln -s /path/to/raw/reads.fastq.gz /path/to/output/working/dir/reads.fastq.gz`


### Step 1: Generate a redundant database of k-mers (RDB)
Generate a k-mer database of all k-mers present in the repeats provided using kmc


Inputs:

1. Repeat annotation bedfile with tab separated columns: [chromosome start_coordinate, stop_coordinate, repeat_type, length, +, repeat_class] where length and + are optional (depending on downstream processing)
2. Genome fasta file
3. Chromosome lengths file with tab separated columns: [chromosome, length]
4. Rootfolder for outputs
5. Maximum memory in GB (ex. 300)

Steps:

1. Create a fasta file for each repeat in annotation file based on genomic coordinates provided
2. Generate a k-mer database for each repeat fasta using KMC
3. Creates a complement fasta including all genomic regions not included in the repeats annotation bedfile using bedtools getfasta and then KMC to generate a background k-mer database.
4. Subtracts background k-mers from each repeat kdb using KMC
5. Dumps background subtracted kdbs and combines all repeat types into one database. This is a redundant k-mer database (RDB) because different repeats can contain the same k-mer (with a different kmer_id). Each fasta ID line is a repeat_type or repeat_class and all of the k-mers found in that repeat type/class.

Outputs:

1. annotations.bed - sorted annotation bedfile from input #1
2. repeats.withClass.txt; repeats.txt - Lists of repeats with/without corresponding class level info
3. Split fasta files (one for each repeat type) & k-mer database
4. repeats.fa - Combined fasta of all repeats
5. annotations.BKG.bed - background coordinates (all regions not in repeat annotation bedfile)
6. Background subtracted k-mer databases for each repeat type
7. Combined redundant k-mer database (RBD.fa)



Run cask_generate_RDB.sh as:

```
cask_generate_RDB.sh {annot_file} {genome_fa} {genome_chromosome_length_file} {kdb_rootfolder} {max mem in gb}

```
Example:
```
cask_generate_RDB.sh "/home/groups/astraigh/kelsey/chm13_genome/v2/censat.clean.chrclass.bed" "/oak/stanford/groups/astraigh/T2T/CURRENT/T2T-CHM13v2.0-hs1_jan2022/hs1.fa" "/home/groups/astraigh/kelsey/chm13_genome/v2/hs1.chromsizes.sorted.txt" "/scratch/groups/astraigh/kelsey/centrochar/cen_only" 350

```


## Step 2: Generate non-redundant k-mer database (NRDB)
Generate non-redundant k-mer database for k-mers only found once in the RDB (unq) and those found more than once (ci2). This gives each k-mer a unique ID that will be used to annotate reads.

Inputs:

1. db type (unq or ci2)
2. rootfolder (same as step 1)
3. Maximum memory allocated in GB (ex. 300)
4. Maximum threads (ex. 24)

Steps:

1. Generates a non-redundant k-mer database (NRDB) from combined repeats.fa where each k-mer receives a unique ID
2. Adds the unique ID from NRDB to the RDB fasta to create an annotated repeat database (ADB)
3. Creates a look up table to map which repeat types or classes each k-mer is found in. Example: kid_2 is found in [repeat_1, repeat_2, repeat_4, repeat_98]

Outputs:

1. NRDB.unq.dump.fa / ci2.dump.fa - fasta file of all 

Run as:

```cask_generate_NRDB.sh {kdb_rootfolder} {unq or ci2} {max mem in gb} {max threads}```

Example:
```
cask_generate_NRDB.sh  unq "/scratch/groups/astraigh/kelsey/centrochar/cen_only" 350 24

cask_generate_NRDB.sh  ci2 "/scratch/groups/astraigh/kelsey/centrochar/cen_only" 350 24
```

## Step 3: Annotate Reads with bbduk using NRDB k-mer IDs
Annotate k-mers in sequencing reads (from fastq) with their unique ID from NRDB

Inputs:

1. k-mer database rootfolder (same as step1 and 2)
2. db type (unq or ci2)
3. Input fastq (can be a symlink)
4. Output subdirectory (directory within the dir where you have raw reads or linked to raw reads in Step 0). No leading `/`
5. Maximum memory allocated in GB (ex. 300)
6. Max threads (ex. 24)

Steps:

1. Creates output directory if it doesn't already exist
2. Runs bbduk using NRDB as reference and fastq as input. Each read name is modified to include the k-mer IDs of NRDB k-mers found in the read as well as their count. 

Outputs:

1. Annotated reads fastq.gz "{outdir_name}/{input_fastq_filename}.ANNOTATED.{dbtype}.fastq.gz" with each read ID line renamed to include the k-mers found in the read as such: 

    ReadID`\t`kmerID_1,kcount_in_NRDB=kcount_in_read` `kmerID_2,kcount_in_NRDB=kcount_in_read

Run as:
 
```cask_annotate_reads.sh {kdb_rootfolder} {db: uq or ci2} {input_fastq} {out_subdir_name} {max_mem} {max_threads}```

Example:
```
cask_annotate_reads.sh "/scratch/groups/astraigh/kelsey/centrochar/cen_only" unq "/scratch/groups/astraigh/kelsey/centrochar/charseq_reads/decon/WT_K562/WTrep1/dna.fastq.gz" "cen_kmatch/k25.minusBKG" 600 24 
```

## Step 4: Process Reads
Use kdb_parser to classify reads into ambivalence groups according to the repeats where the k-mers in the read are found. See CASK_details.md for more specific information.

Inputs:

1. Type list (repeats.withClass.txt output in step 1)
2. kmer LUT from step 2 (LUT.ci2.txt and LUT.unq.txt)
3. Output subdirectory (same as input #4 in step3)
4. dbtype (ci2 or unq)
5. Fastq list file
6. Max number of CPUs

Steps:

1. Processes all reads for each file by comparing k-mer content and the corresponding repeats those k-mers are found in. If only a single repeat type/class is shared amongst all the k-mers, the read is assigned to that repeat type/class which is given a unique, positive ID number. If multiple repeat types/classes are possible, an ambivalence group is created that includes all of those possibilities which is given a unique, negative ID number. Conflicts are given an ID number of -1. If no repeat k-mer are found in the read, the read is given an ID number of 0.

    For example:

    **kmer LUT**


        kmer_ID : [repeats kmer was found in]

        kmer_1 : [repeat_1, repeat_5, repeat_6, repeat_7, repeat_567]

        kmer_5 : [repeat_1, repeat_5, repeat_567]

        kmer_6 : [repeat_1, repeat_5, repeat_23]

        kmer_456 : [repeat_23]

    **Reads:**

        Read 1 contains kmer_1 kmer_5 and kmer_6
        
    Only repeat_1 and repeat_5 are shared across all k-mers, so read is assigned ambiv group ID -2 (the first ambiv group)
        
            


        Read 2 contains kmer_6 and kmer_456
    Only repeat_23 is shared, so read is assigned ambiv group ID 1 (the first non ambiguous assignment)

        Read 3 contains kmer_1, kmer_6

    Possible repeat types are in conflict, read is assigned ambiv group ID -1 

    **Resulting look up dictionary**
    ```
    ambiv_group_ID  : [repeats]
    1 : [repeat_23]
    -2 : [repeat_1, repeat_5]
    ```
        

2. Creates a look up dictionary (as a .pickle file) to map the ambivalence group ID number to the assigned repeat type(s)/class(es)

Outputs:

1. {fastq_filename}.ANNOTATED.{dbtype}.ag.txt
2. {fastq_filename}.ANNOTATED.{dbtype}.aglut.pickle

To run:


Prepare fastq list file (input #5). Create a txt file where each line is the complete path to a fastq you want to process. Note: Output files will be put in the output subdirectory within the base directory from these files. For example, if the fastq_list_file.txt looks like this:
```
/scratch/groups/astraigh/kelsey/centrochar/charseq_reads/decon/ESrep1/dna.fastq.gz
/scratch/groups/astraigh/kelsey/centrochar/charseq_reads/decon/ESrep1/rna.fastq.gz
/scratch/groups/astraigh/kelsey/centrochar/charseq_reads/decon/ESrep2/dna.fastq.gz
/scratch/groups/astraigh/kelsey/centrochar/charseq_reads/decon/ESrep2/rna.fastq.gz
```
Even if these filepaths above are symlinks to another location, outputs will be placed in 
    `/scratch/groups/astraigh/kelsey/centrochar/charseq_reads/decon/ESrep2/subdirectory_from_input_3`

Run as:
```
python /path/to/cask/scripts/cask_process_reads.py {/path/to/kmerdb/repeats.withClass.txt} {/path/to/kmerdb/kdb/k25.minusBKG/LUT.dbtype.txt} {output_subdir} {dbtype: ci2 or unq} {path/to/fastq_list_file.txt} {num_cpus}
```

Example:

``` 
python cask/scripts/cask_process_reads.py /scratch/groups/astraigh/kelsey/centrochar/cen_only/repeats.withClass.txt /scratch/groups/astraigh/kelsey/centrochar/cen_only/kdb/k25.minusBKG/LUT.ci2.txt cen_kmatch/k25.minusBKG ci2 /scratch/groups/astraigh/kelsey/centrochar/charseq_reads/decon/WT_K562/fastq_list.txt 20

python cask/scripts/cask_process_reads.py /scratch/groups/astraigh/kelsey/centrochar/cen_only/repeats.withClass.txt /scratch/groups/astraigh/kelsey/centrochar/cen_only/kdb/k25.minusBKG/LUT.unq.txt cen_kmatch/k25.minusBKG unq /scratch/groups/astraigh/kelsey/centrochar/charseq_reads/decon/WT_K562/fastq_list.txt 20



OR with parallel...

parallel -j2 'python cask/scripts/cask_process_reads.py /scratch/groups/astraigh/kelsey/centrochar/cen_only/backup/cen_only/repeats.withClass.txt /scratch/groups/astraigh/kelsey/centrochar/cen_only/backup/cen_only/kdb/k25.minusBKG/LUT.{}.txt cen_kmatch/k25.minusBKG {} /scratch/groups/astraigh/kelsey/centrochar/charseq_reads/decon/WT_K562/fastq_list.txt 10' ::: ci2 unq

```
## Step 5: Tally annotated reads
Merges processed reads output files for ci2 and unq and counts the number of reads classified in each ambivalence group.

Inputs:

1. Output subdirectory (same as input #4 in step3 and input #3 from step 4)

2. Input fastq file (same as input #3 in step 3).
    
-  Root directory is pulled from this fastq file to get Input files:

     /{fastq file_name}.ANNOTATED.ci2.ag.txt

     {input fastq file directory}/{output subdirectory}/{fastq file_name}.ANNOTATED.unq.ag.txt

Steps:

1. Sets working_directory as {input_fastq_file_directory}/{output_subdirectory}
2. Joins ci2 and unq ag files:
    /{working_directory}/{fastq_file_name}.ANNOTATED.ci2.ag.txt
    /{working_directory}/{fastq_file_name}.ANNOTATED.unq.ag.txt


Outputs:

1. Combined ci2 and unq .ag.txt file (/{working_directory}/{fastq_file_name}.ANNOTATED.ci2_and_unq.ag.txt)

    Output fields:

    -   readID
    - repeat_type_ambivgroup_ci2 (ambiv group assigned based on kmers that appear at least twice in repeat db)
    - class_ambivgroup_ci2
    - number of ci2 kmers in read (kmers that occur at least twice in repeat db)
    - repeat_type_ambivgroup_unq (ambiv group assigned based on kmers that only appear oce in the repeat db)
    - class_ambivgroup_unq 
    - number of unq kmers in read (kmers that occur only once in repeat db)

2. Tallied combined ci2 and unq file (/{working_directory}/{fastq_file_name}.ANNOTATED.ci2_and_unq.tally.txt)

    Output fields:

    - repeat_type_ambivgroup_ci2
    - class_ambivgroup_ci2
    - repeat_type_ambivgroup_unq
    - class_ambivgroup_unq
    - num reads that have this combo of ambiv groups
    - sum of the number of ci2_kmers found in each of those reads
    - sum of the number of unq_kmers found in each of those reads
