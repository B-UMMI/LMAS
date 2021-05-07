#!/bin/bash
# data downloadscript - ZymoBIOMICS Microbial Community Standard

mkdir data

# get description file
wget https://zenodo.org/record/4742651/files/about.md -P data/

# triple reference
mkdir data/reference
wget https://zenodo.org/record/4742651/files/ZymoBIOMICS_genomes.fasta -P data/reference/

# read data
mkdir data/fastq

#  - real even distributed data
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR298/003/ERR2984773/ERR2984773_1.fastq.gz -P data/fastq/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR298/003/ERR2984773/ERR2984773_2.fastq.gz -P data/fastq/

#  - real log distributed data
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR293/005/ERR2935805/ERR2935805_1.fastq.gz -P data/fastq/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR293/005/ERR2935805/ERR2935805_2.fastq.gz -P data/fastq/

#  - mock even distributed data no error
wget https://zenodo.org/record/4742651/files/ENN_1.fq.gz -P data/fastq/
wget https://zenodo.org/record/4742651/files/ENN_2.fq.gz -P data/fastq/

#  - mock even distributed data illumina hiseq error
wget https://zenodo.org/record/4742651/files/EHS_1.fq.gz data/fastq/
wget https://zenodo.org/record/4742651/files/EHS_2.fq.gz data/fastq/

#  - mock log distributed data no error
wget https://zenodo.org/record/4742651/files/LNN_1.fq.gz data/fastq/
wget https://zenodo.org/record/4742651/files/LNN_2.fq.gz data/fastq/

#  - mock log distributed data illumina hiseq error
wget https://zenodo.org/record/4742651/files/LHS_1.fq.gz data/fastq/
wget https://zenodo.org/record/4742651/files/LHS_2.fq.gz data/fastq/