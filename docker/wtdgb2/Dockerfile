# Base image - Alternatives: Conda, Python...
FROM continuumio/miniconda3

LABEL software="wtdbg2" \
    software.version="2.5" \
    about.home="https://github.com/ruanjue/wtdbg2" \
    about.summary="wtdbg2 assembler for long read data" \
    author="Inês Mendes <cimendes@medicina.ulisboa.pt>"

# Maintainer
LABEL maintainer="Inês Mendes <cimendes@medicina.ulisboa.pt>"

# System essentials
RUN conda install -c conda-forge -c bioconda -c defaults wtdbg=2.5

# PS for nextflow
RUN apt-get install -y procps