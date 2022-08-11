FROM continuumio/miniconda:4.7.12

LABEL software="StrainXpress" \
    software.version="11.08.2022" \
    about.home="https://github.com/kangxiongbin/StrainXpress" \
    about.summary="De novo assembler based on overlap-layout-consensus to assemble high complexity metagenome sequencing data at strain resolution." \
    about.documentation="https://github.com/kangxiongbin/StrainXpress#readme" \
    about.license="https://github.com/kangxiongbin/StrainXpress/blob/main/LICENSE" \
    author="Inês Mendes <cimendes@medicina.ulisboa.pt>"

WORKDIR /NGStools

RUN apt-get --allow-releaseinfo-change update \
    && apt-get install procps git build-essential cmake wget libboost-all-dev -y 

RUN conda install -c bioconda python=3.6 scipy pandas minimap2

RUN git clone https://github.com/kangxiongbin/StrainXpress.git \
    && cd StrainXpress \
    && sh install.sh

ENV PATH="/NGStools/StrainXpress/scripts:$PATH"