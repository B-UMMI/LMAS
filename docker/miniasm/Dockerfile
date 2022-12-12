# Base image - Alternatives: Conda, Python...
FROM continuumio/miniconda3

LABEL software="miniasm" \
    software.version="0.3_r179" \
    about.home="https://github.com/lh3/miniasm" \
    about.summary="Ultrafast de novo assembly for long noisy reads (though having no consensus step)" \
    author="Inês Mendes <cimendes@medicina.ulisboa.pt>"

# System essentials
RUN conda install -c bioconda miniasm

# PS for nextflow
RUN apt-get install -y procps