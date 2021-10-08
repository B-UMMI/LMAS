FROM ubuntu:16.04

LABEL software="IDBA" \
    software.version="1.1.3" \
    about.home="https://github.com/loneknightpy/idba" \
    about.summary="Basic iterative de Bruijn graph assembler for second-generation sequencing reads" \
    about.documentation="https://github.com/loneknightpy/idba#readme" \
    author="Inês Mendes <cimendes@medicina.ulisboa.pt>"

RUN  apt-get update && apt-get -y install git \
    gcc \
    autoconf \
    automake \
    g++ \
    make \
    wget


WORKDIR /NGStools

RUN wget https://github.com/loneknightpy/idba/archive/1.1.3.tar.gz
RUN tar xf 1.1.3.tar.gz && rm 1.1.3.tar.gz
# change this parameter so it will allow reads longer than 128 - https://github.com/loneknightpy/idba/issues/26
#RUN sed -i 's#static const uint32_t kMaxShortSequence = 128#static const uint32_t kMaxShortSequence = 256#g' /NGStools/idba-1.1.3/src/sequence/short_sequence.h
#RUN sed -i 's#static const uint32_t kNumUint64 = 4#static const uint32_t kNumUint64 = 8#g' /NGStools/idba-1.1.3/src/basic/kmer.h

RUN cd idba-1.1.3 && sh build.sh

ENV PATH="/NGStools/idba-1.1.3/bin:$PATH"

WORKDIR /data