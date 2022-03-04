# Base image - Alternatives: Conda, Python...
FROM continuumio/miniconda3

LABEL software="ra" \
    software.version="0.2.1" \
    about.home="https://github.com/lbcb-sci/ra" \
    about.summary="Assembler for raw de novo genome assembly of long uncorrected reads." \
    about.license="https://github.com/lbcb-sci/ra/blob/master/LICENSE" \
    author="Inês Mendes <cimendes@medicina.ulisboa.pt>"

# System essentials
RUN conda install -c bioconda ra
RUN apt-get update && apt-get -y install build-essential \
                                         cmake \
                                         zlib1g-dev

# Install RA
RUN git clone --recursive https://github.com/rvaser/ra.git ra
WORKDIR ra
RUN mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release ..
RUN cd build && make

ENV PATH /ra/build/bin:$PATH

# PS for nextflow
RUN apt-get install -y procps