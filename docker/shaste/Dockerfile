# Base image - Alternatives: Conda, Python...
FROM continuumio/miniconda3

LABEL software="shasta" \
    software.version="0.6.0" \
    about.home="https://github.com/chanzuckerberg/shasta" \
    about.summary="De novo assembly from Oxford Nanopore reads." \
    author="Inês Mendes <cimendes@medicina.ulisboa.pt>"

# System essentials
RUN conda install -c bioconda shasta

# PS for nextflow
RUN apt-get install -y procps