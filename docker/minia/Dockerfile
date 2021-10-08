FROM ubuntu:16.04

LABEL software="minia" \
    software.version="3.2.6" \
    about.home="https://github.com/GATB/minia/" \
    about.summary="Minia is a short-read assembler based on a de Bruijn graph." \
    about.documentation="https://github.com/GATB/minia/raw/master/doc/manual.pdf" \
    about.license="https://github.com/GATB/minia/blob/master/LICENSE" \
    author="Inês Mendes <cimendes@medicina.ulisboa.pt>"

RUN apt-get update && apt-get -y install build-essential \
    zlib1g-dev \
    wget \
    python-dev \
    git

WORKDIR /usr/local/src

# Install updated cmake
RUN wget https://github.com/Kitware/CMake/releases/download/v3.18.0/cmake-3.18.0-Linux-x86_64.tar.gz
RUN tar xf cmake-3.18.0-Linux-x86_64.tar.gz && rm cmake-3.18.0-Linux-x86_64.tar.gz
ENV PATH="/usr/local/src/cmake-3.18.0-Linux-x86_64/bin:$PATH"

WORKDIR /NGStools

# Install minia
RUN wget https://github.com/GATB/minia/releases/download/v3.2.6/minia-v3.2.6-bin-Linux.tar.gz
RUN tar xf minia-v3.2.6-bin-Linux.tar.gz && rm minia-v3.2.6-bin-Linux.tar.gz
RUN mv minia-v3.2.6-bin-Linux minia

ENV PATH="/NGStools/minia/bin:$PATH"

WORKDIR /data