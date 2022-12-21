# Base image - Alternatives: Conda, Python...
FROM continuumio/miniconda3

LABEL software="raven-assembler" \
    software.version="1.8.0" \
    about.home="https://github.com/lbcb-sci/raven" \
    about.summary="Raven-assembler for long read data" \
    author="Inês Mendes <cimendes@medicina.ulisboa.pt>"

# System essentials
RUN conda install -c conda-forge -c bioconda -c defaults raven-assembler=1.8.0

# PS for nextflow
RUN apt-get install -y procps