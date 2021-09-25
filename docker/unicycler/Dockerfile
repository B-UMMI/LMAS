FROM ubuntu:16.04

LABEL software="Unicycler" \
    software.version="0.4.9" \
    about.home="https://github.com/rrwick/Unicycler" \
    about.summary="Hybrid assembly pipeline for bacterial genomes." \
    about.documentation="https://github.com/rrwick/Unicycler#readme" \
    about.license="https://github.com/rrwick/Unicycler/blob/main/LICENSE" \
    author="Inês Mendes <cimendes@medicina.ulisboa.pt>"

#install prerequisites including openjdk-11
RUN apt-get update && apt-get install -y software-properties-common &&\
    add-apt-repository -y ppa:openjdk-r/ppa && apt-get update && apt-get install -y \
    python3\
    python3-pip\
    libz-dev\
    bowtie2\
    ncbi-blast+\
    openjdk-11-jre\
    wget\
    bzip2\
    libncurses5-dev\
    libncursesw5-dev\
    libbz2-dev\
    liblzma-dev\
    cmake \
    libcurl4-openssl-dev \
    git

WORKDIR /NGStools

#get unicycler
RUN wget https://github.com/rrwick/Unicycler/archive/refs/tags/v0.4.9.tar.gz &&\
    mkdir unicycler &&\
    tar -xzvf v0.4.9.tar.gz -C unicycler --strip-components 1 &&\
    rm v0.4.9.tar.gz 

#get pilon
RUN wget https://github.com/broadinstitute/pilon/releases/download/v1.24/pilon-1.24.jar &&\
    mkdir pilon &&\
    mv pilon-1.24.jar pilon

#get spades
RUN wget http://cab.spbu.ru/files/release3.15.3/SPAdes-3.15.3-Linux.tar.gz &&\
    mkdir spades &&\
    tar -xzvf SPAdes-3.15.3-Linux.tar.gz -C spades --strip-components 1 &&\
    rm SPAdes-3.15.3-Linux.tar.gz

#get samtools and bcftools
RUN wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 &&\
    wget https://github.com/samtools/bcftools/releases/download/1.13/bcftools-1.13.tar.bz2 &&\
    mkdir samtools &&\
    tar -vxjf samtools-1.13.tar.bz2 -C samtools --strip-components 1 &&\
    rm samtools-1.13.tar.bz2 &&\
    mkdir bcftools &&\
    tar -vxjf bcftools-1.13.tar.bz2 -C bcftools --strip-components 1 &&\
    rm bcftools-1.13.tar.bz2 &&\
    cd samtools &&\
    make &&\
    make install &&\
    cd ../bcftools &&\
    make &&\
    make install

#get racon
RUN wget https://github.com/isovic/racon/releases/download/1.4.3/racon-v1.4.3.tar.gz &&\
    mkdir racon &&\
    tar -xzvf racon-v1.4.3.tar.gz -C racon --strip-components 1 &&\
    rm racon-v1.4.3.tar.gz &&\
    cd racon &&\
    mkdir build &&\
    cd build &&\
    cmake -DCMAKE_BUILD_TYPE=Release ../ &&\
    make &&\
    make install

#get miniasm
RUN wget https://github.com/lh3/miniasm/archive/v0.3.tar.gz &&\
    mkdir miniasm &&\
    tar -xzvf v0.3.tar.gz -C miniasm --strip-components 1 &&\
    rm v0.3.tar.gz &&\
    cd miniasm &&\
    make

#set env path variable for installed programs
ENV PATH="/NGStools/pilon:${PATH}"
ENV PATH="/NGStools/spades/bin:${PATH}"
ENV PATH="/NGStools/miniasm:${PATH}"

#build unicycler
RUN python3 /NGStools/unicycler/setup.py install

WORKDIR /data