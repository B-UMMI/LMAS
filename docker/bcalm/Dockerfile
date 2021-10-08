FROM ubuntu:16.04

LABEL software="BCALM" \
    software.version="2.2.3" \
    about.home="https://github.com/GATB/bcalm" \
    about.summary="Compacted de Bruijn graph construction in low memory" \
    about.documentation="https://github.com/GATB/bcalm#readme" \
    about.license="https://github.com/GATB/bcalm/blob/master/LICENSE" \
    author="Inês Mendes <cimendes@medicina.ulisboa.pt>"

RUN apt-get update && apt-get -y install build-essential \
    wget \
    git

RUN apt-get -y install software-properties-common
RUN add-apt-repository ppa:ubuntu-toolchain-r/test
RUN apt-get update
RUN apt-get -y install gcc-4.9
RUN apt-get -y upgrade libstdc++6
WORKDIR /NGStools

# Dependencies
# BCALM
RUN wget https://github.com/GATB/bcalm/releases/download/v2.2.3/bcalm-binaries-v2.2.3-Linux.tar.gz
RUN tar xf bcalm-binaries-v2.2.3-Linux.tar.gz && rm bcalm-binaries-v2.2.3-Linux.tar.gz

ENV PATH="/NGStools/bcalm-binaries-v2.2.3-Linux/bin:${PATH}"

WORKDIR /data