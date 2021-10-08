FROM cimendes/lmas-python-base:0.1

LABEL software="minimap2" \
    software.version="2.22" \
    about.home="https://github.com/lh3/minimap2" \
    about.summary="A versatile pairwise aligner for genomic and spliced nucleotide sequences" \
    about.documentation="https://lh3.github.io/minimap2/minimap2.html" \
    about.license="https://github.com/lh3/minimap2/blob/master/LICENSE.txt" \
    author="Inês Mendes <cimendes@medicina.ulisboa.pt>"

RUN apt-get update
RUN apt-get -y install \
    wget \
    lbzip2

WORKDIR /NGStools

RUN wget https://github.com/lh3/minimap2/releases/download/v2.22/minimap2-2.22_x64-linux.tar.bz2 && \
    tar xf minimap2-2.22_x64-linux.tar.bz2 && rm minimap2-2.22_x64-linux.tar.bz2

ENV PATH="/NGStools/minimap2-2.22_x64-linux:$PATH"