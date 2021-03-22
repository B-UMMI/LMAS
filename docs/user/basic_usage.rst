Basic Usage
===========

After you have a local installation of LMAS, the mock community data needs to be downloaded.

The triple reference sequences can be passed with the ``--reference`` parameter, and ``--fastq`` recieves 
the raw data for assembly. The raw data is a collection of sequence fragments from the references, and can 
be either obtained *in silico* or from real sequencing platforms.

The optional parameter ``--md`` allows the user to pass information on input samples to be presented in the report. 


Download ZymoBIOMICS Microbial Community Standard Data
------------------------------------------------------

As proof-of-concept,the eight bacterial genomes and four plasmids of the 
`ZymoBIOMICS Microbial Community Standards <https://www.zymoresearch.com/collections/zymobiomics-microbial-community-standards>`_ 
were used as reference. Raw sequence data of the mock communities, with an even and logarithmic distribution of species, 
and a simulated sample of the evenly distributed reads generated from the genomes. 

The triple-reference sequences and the mock sample are available at zenodo: https://zenodo.org/record/4588970#.YEeA83X7RhE

The even and log distributed raw sequence data is available at https://www.ebi.ac.uk/ena/browser/view/ERR2984773 and 
https://www.ebi.ac.uk/ena/browser/view/ERR2935805, respectively. 

A script to download and structure the ZymoBIOMICS data to be ready to used as default input for LMAS is provided, 
included in LMAS repository. To run it, simply execute: 

.. code-block:: bash

    sh get_data.sh 

The files will be saved in the following structure: 

.. code-block:: bash

    data/
    ├── fastq
    │   ├── ERR2935805_1.fq.gz 
    │   ├── ERR2935805_2.fq.gz
    │   ├── ERR2984773_1.fq.gz
    │   ├── ERR2984773_2.fq.gz
    │   ├── mockSample_1.fq.gz
    │   └── mockSample_2.fq.gz
    └── reference
        └── Zymos_Genomes_triple_chromosomes.fasta


Customizing LMAS
----------------

Users can customize the workflow execution either by using command line options or by modifying a simple 
plain-text configuration file (``params.config``), where parameters are set as key-value pairs. 

The version of tools used can also be changed by providing new container tags in the appropriate configuration file 
(``containers.config``), as well as the resources for each process (``resources.config``).

