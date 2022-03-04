# Base image - Alternatives: Conda, Python...
FROM continuumio/miniconda3

LABEL software="Flye" \
    software.version="2.9" \
    about.home="https://github.com/fenderglass/Flye" \
    about.summary="Flye Assembler for long read data" \
    author="Inês Mendes <cimendes@medicina.ulisboa.pt>"

# System essentials
RUN conda install -c conda-forge -c bioconda -c defaults flye=2.9

# Set up workdir for installation 
WORKDIR /data

# PS for nextflow
RUN apt-get install -y procps