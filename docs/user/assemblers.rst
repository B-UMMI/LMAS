Short-Read (Meta)Genomic Assemblers
===================================

We've compiled a collection of *de novo* assembly tools, including **Overlap, Layout and Consensus (OLC)** 
and **De Bruijn graph** assembly algorithms, both **single k-mer and multiple k-mer value approaches**, and hybrid assemblers.
The collection includes both genomic and metagenomic assemblers developed explicitly to handle metagenomic datasets.

Selection Criteria
-------------------

Open-source tools, with clear documentation describing the methodology implemented, were favoured. 
The collection of tools were ordered by the date of the last update and a Docker container for the top 12 assemblers 
was created with the latest released version, with the version used as the tag. 
In the case of tools with no release, the container was compiled with the latest version in the default branch of the 
source repository, using the date of the last update as the tag.

Assemblers in LMAS
------------------

In its current form, 12 assemblers are implemented in LMAS.
To change the version of a particular assembler is as simple as substituting the container for the process 
of that assembler in LMAS’ container configuration file (see https://lmas.readthedocs.io/en/latest/user/basic_usage.html#containers-config).

Assemblers benchmarked in LMAS, in alphabetical order:

AbYSS
^^^^^^

This assembler, published by `Jackman et al, 2017 <http://doi.org/10.1101/gr.214346.116>`_ in 
*Genome Research*, is a *de novo* sequence assembler intended for short paired-end reads and genomes of all sizes.
It follows the model of Minia, wherein a probabilistic Bloom filter representation is used to encode the de Bruijn 
graph, reducing memory requirements for de novo assembly. **It's a traditional single k-mer value De Bruijn assembler.**

* **Source code:** https://github.com/bcgsc/abyss
* **Date of last release:** 22/04/2021
* **Container:** `cimendes/abyss:2.3.1-1 <https://hub.docker.com/repository/docker/cimendes/abyss>`_ 

GATB-Minia Pipeline
^^^^^^^^^^^^^^^^^^^

GATB-Minia is an assembly pipeline, still unpublished, that consists of Bloocoo for error correction, Minia 3 for contigs 
assembly, which is based on the BCALM2 assembler, and BESST for scaffolding.
It was developed to extend the Minia assembler to use **De Bruijn algorithm with multiple** **k-mer values**.
**It was developed explicitly to handle metagenomic data.**

* **Source code:** https://github.com/GATB/gatb-minia-pipeline
* **Date of last release:** 31/07/2020
* **Container:** `cimendes/gatb-minia-pipeline:31.07.2020-1 <https://hub.docker.com/repository/docker/cimendes/gatb-minia-pipeline>`_

IDBA-UD
^^^^^^^ 

Published by `Peng et al. 2012 <https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/bts174>`_, it's 
a **De Bruijn graph assembler for assembling reads from single-cell sequencing or metagenomic sequencing technologies** with 
uneven sequencing depths. It employs multiple depth relative thresholds to remove erroneous k-mers in both low-depth and 
high-depth regions. The technique of local assembly with paired-end information is used to solve the branch problem of 
low-depth short repeat regions. To speed up the process, an error correction step is conducted to correct reads of 
high-depth regions that can be aligned to high confidence contigs.

* **Source code:** https://github.com/loneknightpy/idba
* **Date of last release:** 11/07/2016
* **Container:** `cimendes/idba:1.1.3-1 <https://hub.docker.com/repository/docker/cimendes/idba>`_

MEGAHIT
^^^^^^^

MEGAHIT, published by `Li et al. 2015 <https://academic.oup.com/bioinformatics/article/31/10/1674/177884>`_, is a 
*de novo* assembler for **assembling large and complex metagenomics data** in a time- and cost-efficient manner. 
It makes use of the succinct **de Bruijn graph, with a multiple k-mer size** strategy. In each iteration, MEGAHIT cleans 
potentially erroneous edges by removing tips, merging bubbles and removing low local coverage edges, especially 
useful for metagenomics which suffers from non-uniform sequencing depths.

* **Source code:** https://github.com/voutcn/megahit
* **Date of last release:** 15/10/2019
* **Container:** `cimendes/megahit-assembler:1.2.9-1 <https://hub.docker.com/repository/docker/cimendes/megahit-assembler>`_

MetaHipMer2
^^^^^^^^^^^^

MetaHipMer2, published by `Georganas et al. 2018 <https://doi.org/10.1109/SC.2018.00013>`_, is a 
a high-quality and high-performance metagenome assembler that employs an iterative de Bruijn graph approach. 
MetaHipMer leverages a specialized scaffolding algorithm that produces long scaffolds and accommodates 
the idiosyncrasies of metagenomes.

* **Source code:** https://bitbucket.org/berkeleylab/mhm2
* **Date of last release:** 13/10/2020
* **Container:** `cimendes/mhm2:v2.0.0-65-gaad446d-generic <https://hub.docker.com/repository/docker/cimendes/mhm2>`_

METASPADES
^^^^^^^^^^

SPAdes started as a tool aiming to resolve uneven coverage in single-cell genome data, but later metaSPAdes 
was released, building a **specific metagenomic pipeline on top of SPAdes**. It was published by `Nurk et al. 2017 <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5411777/>`_, 
and like SPAdes, it uses **multiple k-mer sizes of the De Bruijn graph**, starting with the lowest k-mer size and adding 
hypothetical k-mers to connect the graph.

* **Source code:** https://github.com/ablab/spades
* **Date of last release:** 23/07/2021
* **Container:** `cimendes/spades:3.15.3-1 <https://hub.docker.com/repository/docker/cimendes/spades>`_


MINIA
^^^^^

This tool, published by `Chikhi & Rizk, 2013 <https://almob.biomedcentral.com/articles/10.1186/1748-7188-8-22>`_ in 
*Algorithms for Molecular Biology*, performs the assembly on a data structure based on unitigs produced by the BCALM 
software and using graph simplifications that are heavily inspired by the SPAdes assembler. Minia is a short-read 
traditional assembler based on **De Bruijn graph using a single k-mer length**.

* **Source code:** https://github.com/GATB/minia
* **Date of last release:** 26/05/2021
* **Container:** `cimendes/minia:3.2.6-1 <https://hub.docker.com/repository/docker/cimendes/minia>`_ 

SKESA
^^^^^

This *de novo* sequence read assembler is based on **De Bruijn graphs** and uses conservative heuristics and is designed 
to create breaks at repeat regions in the genome, creating shorter assemblies but with greater sequence quality. It 
tries to obtain good contiguity by using **multiple k-mers** longer than mate length and up to insert size. It was published by 
`Souvorov et al. 2018 <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1540-z>`_. 

* **Source code:** https://github.com/ncbi/SKESA
* **Date of last update:** 02/04/2021
* **Container:** `cimendes/skesa:2.5.0-1 <https://hub.docker.com/repository/docker/cimendes/skesa>`_

SPADES
^^^^^^

A tool aiming to resolve uneven coverage in **single-cell genome data through multiple k-mer sizes of De Bruijn graphs**. 
It starts with the smallest k-mer size and adds hypothetical k-mers to connect the graph. It was published by
`Bankevich et al. 2012 <https://pubmed.ncbi.nlm.nih.gov/22506599/>`_. 

* **Source code:** https://github.com/ablab/spades
* **Date of last release:** 23/07/2021
* **Container:** `cimendes/spades:3.15.3-1 <https://hub.docker.com/repository/docker/cimendes/spades>`_

UNICYCLER
^^^^^^^^^

An assembly pipeline for **bacterial genomes** that can do long-read assembly, hybrid assembly and short-read assembly. 
When assembling Illumina-only read sets where it functions as a SPAdes-optimiser, using a **de Bruijn algorithm with** 
**multiple k-mer values**. It was published by `Wick et al. 2017 <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005595>`_.

* **Source code:** https://github.com/rrwick/Unicycler
* **Date of last release:** 03/05/2021
* **Container:** `cimendes/unicycler:0.4.9-1 <https://hub.docker.com/repository/docker/cimendes/unicycler>`_

VELVETOPTIMISER
^^^^^^^^^^^^^^^

This optimizing pipeline, developed by Torsten Seeman, is still unpublished but extends the original Velvet assembler by 
performing **several de Bruijn assemblies with variable k-mer sizes**. It searches a supplied hash value range for the optimum, estimates 
the expected coverage and then searches for the optimum coverage cutoff. It uses Velvet's internal mechanism for estimating 
insert lengths for paired-end libraries. It can optimise the assemblies by either the default optimisation condition or by a 
user-supplied one. It outputs the results to a subdirectory and records all its operations in a logfile.

* **Source code:** https://github.com/tseemann/VelvetOptimiser
* **Date of last update:** 21/01/2017
* **Container:** `cimendes/velvetoptimiser:2.2.6-1 <https://hub.docker.com/repository/docker/cimendes/velvetoptimiser>`_
