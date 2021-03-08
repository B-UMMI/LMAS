#!/bin/bash
# data downloadscript - ZymoBIOMICS Microbial Community Standard

mkdir data

# triple reference
mkdir data/reference
wget https://zenodo.org/record/4588970/files/Zymos_Genomes_triple_chromosomes.fasta -P data/reference/

# read data
mkdir data/fastq
#  - mock data
wget https://zenodo.org/record/4588970/files/mockSample_1.fq.gz -P data/fastq/
wget https://zenodo.org/record/4588970/files/mockSample_2.fq.gz -P data/fastq/

#  - even distributed data
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR298/003/ERR2984773/ERR2984773_1.fastq.gz -P data/fastq/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR298/003/ERR2984773/ERR2984773_2.fastq.gz -P data/fastq/

#  - log distributed data
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR293/005/ERR2935805/ERR2935805_1.fastq.gz -P data/fastq/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR293/005/ERR2935805/ERR2935805_2.fastq.gz -P data/fastq/