FROM ubuntu:16.04

LABEL software="SPAdes" \
    software.version="3.15.3" \
    about.home="https://github.com/ablab/spades" \
    about.summary="SPAdes Genome Assembler" \
    about.documentation="https://github.com/ablab/spades#readme" \
    author="Inês Mendes <cimendes@medicina.ulisboa.pt>"

RUN  apt-get update && apt-get -y install git \
    wget \
    python \
    build-essential \
    cmake

RUN apt-get -y install zlib1g-dev libbz2-dev

WORKDIR /NGStools

RUN wget https://github.com/ablab/spades/releases/download/v3.15.3/SPAdes-3.15.3-Linux.tar.gz
RUN tar xf SPAdes-3.15.3-Linux.tar.gz && rm SPAdes-3.15.3-Linux.tar.gz
RUN mv SPAdes-3.15.3-Linux spades
ENV PATH="/NGStools/spades/bin:$PATH"