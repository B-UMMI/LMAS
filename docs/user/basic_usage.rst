Basic Usage
===========

To use LMAS, the **reference sequences** must be passed with the ``--reference`` parameter, and ``--fastq`` 
receives the **short-read paired-end raw data** for assembly. 

The optional parameter ``--md`` allows the user to pass information on input samples, in a markdown file, to be 
presented in the LMAS report. 

All complete genomes (reference linear replicons) should be provided in a single file. 
Due to the ambiguous starting position of a circular replicon, an assembled contig will typically not align to 
the reference in a single unbroken alignment. Therefore, the linearized replicons are concatenated 
three times by LMAS to ensure that contigs can fully align even with start-end overlap and regardless of 
their starting position relative to that of the reference. 

The raw data is a collection of sequence fragments 
from the references and can be either obtained *in silico* or from real sequencing platforms.

.. warning:: By default, LMAS expects the input data in a ``data/`` folder, with the reference sequences in ``data/reference/*.fasta``, the read data in ``data/fastq/*_{1,2}.*``, and the markdown file in ``data/*.md``.

When you clone it, LMAS has the following folder structure:

.. code-block:: bash

    LMAS
    ├── bin/
    ├── containers.config
    ├── docker/
    ├── docs/
    ├── get_data.sh
    ├── lib/
    ├── LICENSE
    ├── LMAS.nf
    ├── nextflow.config
    ├── params.config
    ├── profiles.config
    ├── README.md
    ├── resources/
    ├── resources.config
    └── templates/

* The ``LMAS.nf`` is the main execution file for LMAS. 
* The ``get_data.sh`` bash script file downloads the ZymoBIOMICS Microbial Community Standard data.
* The ``containers.config``, ``nextflow.config``, ``params.config``, ``profiles.config`` and ``resources.config`` are LMAS configuration files.
* The ``bin/`` and ``templates/`` folders contain custom LMAS code for data processing.
* The ``docs/`` folder contains LMAS documentation source files.
* The ``docker/`` folder contains the dockerfile for LMAS' base container.
* The ``resources/`` folder contains the LMAS report compiled code.


Customizing LMAS workflow configuration
---------------------------------------

Users can customize the **workflow execution** either by using **command-line options**, with ``--<name of parameter> <option>``
or by modifying a simple **plain-text configuration file**, where parameters are set as key-value pairs.

There are four configuration files in LMAS:

nextflow.config
^^^^^^^^^^^^^^^

This is Nextflow main configuration file. **It should not be edited**. 

params.config
^^^^^^^^^^^^^

The ``params.config`` file includes all available parameters for LMAS and their respective default values.

containers.config 
^^^^^^^^^^^^^^^^^

The ``containers.config`` file includes the container directive for each process in LMAS. 
These containers are retrieved from **dockerhub** if they do not exist locally yet. 

.. warning:: You can change the container string to any other value, but it should point to an image that exists on dockerhub or locally.

profiles.config 
^^^^^^^^^^^^^^^

The ``profiles.config`` file includes a set of pre-made profiles with all possible combinations of executors and container engines. 
You can add new ones or modify an existing one.

resources.config 
^^^^^^^^^^^^^^^^
 
The ``resources.config`` file includes the **CPUs** and **memory** directives provided for each assembler in LMAS. 

.. warning:: The **memory** directive increments automatically when a task is retried. If the directive is set to ``{16.Gb*task.attempt}``, the memory used will be 16 Gb multiplied by the number of attempts. 


ZymoBIOMICS Microbial Community Standard Data
-------------------------------------------------

As a proof-of-concept, the eight bacterial genomes and four plasmids of the 
`ZymoBIOMICS Microbial Community Standards <https://www.zymoresearch.com/collections/zymobiomics-microbial-community-standards>`_ 
were used as reference. Raw sequence data of the mock communities, with an even and logarithmic distribution of species, 
and a simulated sample of the evenly distributed reads were generated from the genomes. 

The triple-reference sequences and the mock sample are available at zenodo: https://zenodo.org/record/4588970#.YEeA83X7RhE

The even and log distributed raw sequence data is available at https://www.ebi.ac.uk/ena/browser/view/ERR2984773 and 
https://www.ebi.ac.uk/ena/browser/view/ERR2935805, respectively. 

A script to download and structure the ZymoBIOMICS data to be used as **default input** for LMAS is provided, 
included in LMAS' repository. To run it, simply execute: 

.. code-block:: bash

    sh get_data.sh 

The files will be saved in the following structure: 

.. code-block:: bash

    data/
    ├── about.md
    ├── fastq
    │   ├── ERR2935805_1.fq.gz
    │   ├── ERR2935805_2.fq.gz
    │   ├── ERR2984773_1.fq.gz
    │   ├── ERR2984773_2.fq.gz
    │   ├── mockSample_1.fq.gz
    │   └── mockSample_2.fq.gz
    └── reference
        └── Zymos_Genomes_triple_chromosomes.fasta
        
This is already the expected input for LMAS. To execute LMAS you simply need to call the ``LMAS.nf`` execution file with Nextflow.

.. code-block:: bash

    nextflow run LMAS.nf
