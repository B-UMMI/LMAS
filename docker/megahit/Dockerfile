FROM ubuntu:16.04

LABEL software="MEGAHIT" \
    software.version="1.2.9" \
    about.home="https://github.com/voutcn/megahit" \
    about.summary="Ultra-fast and memory-efficient (meta-)genome assembler" \
    about.documentation="https://github.com/voutcn/megahit#readme" \
    about.license="https://github.com/voutcn/megahit/blob/master/LICENSE" \
    author="Inês Mendes <cimendes@medicina.ulisboa.pt>"
RUN apt-get update

WORKDIR /NGStools

# -- General Dependencies ---
RUN apt-get install -y wget python python3 zlib1g-dev build-essential

# -- megahit ---
RUN wget https://github.com/voutcn/megahit/releases/download/v1.2.9/MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz
RUN tar -xzf MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz && rm MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz
RUN mv MEGAHIT-1.2.9-Linux-x86_64-static megahit

ENV PATH="/NGStools/megahit:$PATH"