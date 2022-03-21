FROM continuumio/miniconda3

LABEL software="CANU" \
    software.version="2.2" \
    about.home="https://github.com/marbl/canu" \
    about.summary="A single molecule sequence assembler for genomes large and small." \
    about.documentation="https://canu.readthedocs.io/en/latest/index.html" \
    about.license="https://github.com/marbl/canu/blob/master/README.licenses" \
    author="Inês Mendes <cimendes@medicina.ulisboa.pt>"

# System essentials
RUN conda install -c conda-forge -c bioconda -c defaults canu=2.2

# Set up workdir for installation 
WORKDIR /data

# PS for nextflow
RUN apt-get install -y procps