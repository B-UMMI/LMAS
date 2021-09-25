FROM ubuntu:16.04

LABEL software="skesa" \
    software.version="2.5.0" \
    about.home="https://github.com/ncbi/SKESA" \
    about.summary="SKESA assembler" \
    about.documentation="https://github.com/ncbi/SKESA#readme" \
    about.license="https://github.com/ncbi/SKESA/blob/master/LICENSE" \
    author="Inês Mendes <cimendes@medicina.ulisboa.pt>"

RUN apt-get update && apt-get -y install \
	bzip2 \
	git \
	python3 \
	wget

RUN apt-get -y install software-properties-common
RUN add-apt-repository ppa:ubuntu-toolchain-r/test
RUN apt-get update
RUN apt-get -y install gcc-4.9
RUN apt-get -y upgrade libstdc++6

WORKDIR /NGStools

RUN wget https://github.com/ncbi/SKESA/releases/download/skesa.2.4.0_saute.1.3.0_1/skesa.centos.7.9.2009
RUN mv skesa.centos.7.9.2009 skesa && chmod +x skesa

ENV PATH="/NGStools:$PATH"