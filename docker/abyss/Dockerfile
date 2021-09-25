FROM continuumio/miniconda:4.7.12

LABEL software="abyss" \
    software.version="2.3.1" \
    about.home="https://github.com/bcgsc/abyss" \
    about.summary="Assemble large genomes using short reads" \
    about.documentation="https://github.com/bcgsc/abyss#readme" \
    about.license="https://github.com/bcgsc/abyss/blob/master/LICENSE" \
    author="Inês Mendes <cimendes@medicina.ulisboa.pt>"

WORKDIR /NGStools

RUN conda install -c bioconda abyss
