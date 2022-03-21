FROM ubuntu:16.04

LABEL software="CANU" \
    software.version="2.2" \
    about.home="https://github.com/marbl/canu" \
    about.summary="A single molecule sequence assembler for genomes large and small." \
    about.documentation="https://canu.readthedocs.io/en/latest/index.html" \
    about.license="https://github.com/marbl/canu/blob/master/README.licenses" \
    author="Inês Mendes <cimendes@medicina.ulisboa.pt>"

RUN apt-get update

WORKDIR /NGStools

# -- General Dependencies ---
RUN apt-get install -y wget curl python python3 zlib1g-dev build-essential procps

# -- canu ---
RUN curl -L https://github.com/marbl/canu/releases/download/v2.2/canu-2.2.Linux-amd64.tar.xz --output canu-2.2.Linux.tar.xz &&\
    tar -xJf canu-2.2.*.tar.xz

ENV PATH="/NGStools/canu-2.2/bin:$PATH"