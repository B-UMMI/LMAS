FROM ubuntu:16.04

LABEL software="PANDAseq" \
    software.version="2.11" \
    about.home="https://github.com/neufeld/pandaseq" \
    about.summary="PAired-eND Assembler for DNA sequences" \
    about.documentation="https://github.com/neufeld/pandaseq#readme" \
    about.license="https://github.com/neufeld/pandaseq/blob/master/COPYING" \
    author="Inês Mendes <cimendes@medicina.ulisboa.pt>"

RUN apt-get update && apt-get -y install build-essential \
    libtool automake zlib1g-dev libbz2-dev pkg-config git

WORKDIR /NGStools

RUN git clone https://github.com/neufeld/pandaseq.git
RUN cd pandaseq && ./autogen.sh && ./configure && make && make install

ENV PATH="/NGStools/pandaseq:$PATH"