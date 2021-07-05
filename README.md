# CASK
Classification of Ambivalent Sequences using K-mers.

CASK is a tool to classify sequencing reads from repetitive genomic regions based on their k-mer composition. CASK works by first generating an annotated k-mer database, linking all the k-mers in the genome with an ambivalence group, which is the set of all repeat types in which the k-mer occurs. CASK then uses the annotated k-mer database to parse sequencing reads and annotate them with an ambivalence group consistent with all the kmers in the read. 

CASK relies on KMC for generating repeat-specific k-mer databases and on Bbduk for extracting relevant kmers from the reads. Generation of the annotated k-mer databases and annotation of the reads is currently done through bash scripts and a small python module.

The inputs for CASK are :
- a fasta file with the genome sequence
- a bed file with the genomic coordinates of the repeats and their type (obtained typically with RepeatMasker)
- a fastq file with sequencing reads to be classified

## Dependencies :
- bedtools (https://bedtools.readthedocs.io/en/latest/)
- KMC (https://github.com/refresh-bio/KMC)
- bbduk.sh from the BBMap suite (https://sourceforge.net/projects/bbmap/) 
- Standard Unix tools
