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
    ├── CITATION.cff
    ├── conf
    ├── docker
    ├── docs
    ├── get_data.sh
    ├── lib
    ├── LICENSE
    ├── LMAS
    ├── main.nf
    ├── modules
    ├── nextflow.config
    ├── README.md
    ├── resources
    ├── templates
    └── test

* The ``LMAS`` is the main execution file for LMAS.
* The ``main.nf`` is the workflow execution file for `Nextflow <https://www.nextflow.io/>`_.
* The ``modules`` folder contains the LMAS `Nextflow <https://www.nextflow.io/>`_  DSL2 modules.
* The ``nextflow.config`` and the files in ``conf/`` are the LMAS configuration files.
* The ``lib/`` and ``templates/`` folders contain custom LMAS code for data processing.
* The ``docs/`` folder contains the LMAS documentation source files.
* The ``docker/`` folder contains the dockerfile for the base container and all assemblers in LMAS.
* The ``resources/`` folder contains the LMAS report compiled code.
* The ``test`` folder contains LMAS test data and files.
* The ``get_data.sh`` bash script file downloads the ZymoBIOMICS Microbial Community Standard data.


Customizing LMAS workflow configuration
---------------------------------------

Users can customize the **workflow execution** either by using **command-line options**, with ``--<name of parameter> <option>``
or by modifying a simple **plain-text configuration file**, located in ``conf/`` folder, where parameters are set as key-value pairs.

We advise to use the command line options directly, and more information is available the in `Parameters <../user/parameters.html>`_ docummentation page.

Unsing the command line options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To use LMAS the following options are available:

.. code-block:: bash

                      _    __  __   _   ___
       /\︵︵/\      | |  |  \/  | /_\ / __|
      (◕('人')◕)     | |__| |\/| |/ _ \\__ \
         |︶|        |____|_|  |_/_/ \_\___/

         Last Metagenomic Assembler Standing

      Input parameters:
         --fastq                    Path expression to paired-end fastq files.
                                    (default: data/fastq/*_{1,2}.*)
         --reference                Path to the genome reference fasta file.
                                    (default: data/reference/*.fasta)
         --md                       Path to markdown with input sample description for report (optional).
                                    (default: data/*.md)

      Mapping and filtering paramenters:
         --minLength                Value for minimum contig length, in basepairs.
                                    (default: 1000)
         --mapped_reads_threshold   Value for the minimum percentage of a read aligning to the
                                    contig to be considered as mapped.
                                    (default: 0.75)

      Assembly quality assessment parameters:
         --n_target                 Target value for the N, NA and NG metrics, ranging from 0 to 1.
                                    (default: 0.5)
         --l_target                 Target value for the L metric, ranging from 0 to 1.
                                    (default: 0.5)
         --plot_scale               Scale of x-axis for the L, NA and NG metrics plots.
                                    Allowed values: 'linear' or 'log'.
                                    (default: log)

      Assembly execution parameters:
         --abyss                    Boolean controling the execution of the ABySS assembler.
                                    (default: true)
         --abyssKmerSize            K-mer size for the ABySS assembler, as an intiger.
                                    (default 96)
         --abyssBloomSize           Bloom filter size for the ABySS assembler.
                                    It must be a sting with a value and an unit.
                                    (default: 2G)
         --bcalm                    Boolean controling the execution of the BCALM2 assembler.
                                    (default: true)
         --bcalmKmerSize            K-mer size for the BCALM2 assembler, as an intiger.
                                    (default 31)
         --gatb_minia               Boolean controling the execution of the GATB Minia Pipeline assembler.
                                    (default: true)
         --gatbKmerSize             K-mer sizes for the GATB Minia Pipeline assembler.
                                    It must be a sting with the values separated with a comma.
                                    (default 21,61,101,141,181)
         --gatb_besst_iter          Number of iteration during Besst scaffolding for the
                                    GATB Minia Pipeline assembler.
                                    (default 10000)
         --gatb_error_correction    Boolean to control weather to skip error correction for the
                                    GATB Minia Pipeline assembler.
                                    (default false)
         --idba                     Boolean controling the execution of the IDBA-UD assembler.
                                    (default true)
         --metahipmer2              Boolean controling the execution of the MetaHipMer2 assembler.
                                    (default true)
         --metahipmer2KmerSize      K-mer sizes for the MetaHipMer2 assembler.
                                    It must be a sting with the values separated with a comma.
                                    (default 21,33,55,77,99)
         --minia                    Boolean controling the execution of the minia assembler.
                                    (default: true)
         --miniaKmerSize            K-mer size for the minia assembler, as an intiger.
                                    (default 31)
         --megahit                  Boolean controling the execution of the MEGAHIT assembler.
                                    (default true)
         --megahitKmerSize          K-mer sizes for the MEGAHIT assembler.
                                    It must be a sting with the values separated with a comma.
                                    (default 21,29,39,59,79,99,119,141)
         --metaspades               Boolean controling the execution of the metaSPAdes assembler.
                                    (default true)
         --metaspadesKmerSize       K-mer sizes for the metaSPAdes assembler.
                                    It must be a sting with 'auto' or the values separated with a space.
                                    (default auto)
         --spades                   Boolean controling the execution of the SPAdes assembler.
                                    (default true)
         --spadesKmerSize           K-mer sizes for the SPAdes assembler.
                                    It must be a sting with 'auto' or the values separated with a space.
                                    (default auto)
         --skesa                    Boolean controling the execution of the SKESA assembler.
                                    (default true)
         --unicycler                Boolean controling the execution of the Unicycler assembler.
                                    (default true)
         --velvetoptimiser          Boolean controling the execution of the VelvetOptimiser assembler.
                                    (default: true)
         --velvetoptimiser_hashs    Starting K-mer size for the VelvetOptimiser assembler, as an intiger.
                                    (default 19)
         --velvetoptimiser_hashe    End K-mer size for the VelvetOptimiser assembler, as an intiger.
                                    (default 31)

      Execution resources parameters:
         --cpus                     Number of CPUs for the assembly and mapping processes, as an intiger.
                                    This resource is double for each retry until max_cpus is reached.
                                    (default 8)
         --memory                   Memory for the assembly and mapping processes, in the format of
                                    'value'.'unit'.
                                    This resource is double for each retry until max_memory is reached.
                                    (default 32 GB)
         --time                     Time limit for the assembly and mapping processes, in the format of
                                    'value'.'unit'.
                                    This resource is double for each retry until max_time is reached.
                                    (default 1d)
         --max_cpus                 Maximum number of CPUs for the assembly and mapping processes,
                                    as an intiger. It overwrites the --cpu parameter.
                                    (default 32)
         --max_memory               Maximum memory for the assembly and mapping processes, in the format of
                                    'value'.'unit'. It overwrites the --memory parameter.
                                    (default 100 GB)
         --max_time                 Maximum time for the assembly and mapping processes, in the format of
                                    'value'.'unit'. It overwrites the --time parameter.
                                    (default 3d)


LMAS configuration files
^^^^^^^^^^^^^^^^^^^^^^^^^

There are four configuration files in LMAS:

nextflow.config
^^^^^^^^^^^^^^^

This is Nextflow main configuration file.
The resource parameters are available here, and can be changed directly here. 

.. warning:: The **memory** and **cpu** directives increment automatically when a task is retried. If the directive is set to ``{16.Gb*task.attempt}``, the memory used will be 16 Gb multiplied by the number of attempts. By default LMAS is set to run a maximum of 2 retires per process. If the maximum resources are reached before the maximum number of tries, these won't be incremented beyond the defined limit.

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


ZymoBIOMICS Microbial Community Standard Data
-------------------------------------------------

As a proof-of-concept, the eight bacterial genomes and four plasmids of the 
`ZymoBIOMICS Microbial Community Standards <https://www.zymoresearch.com/collections/zymobiomics-microbial-community-standards>`_ 
were used as reference. Raw sequence data of the mock communities, with an even and logarithmic distribution of species both from a 
real sequencing run and a simulated dataset, with and without error, matching the real data distribution of species, were used as input for LMAS. 

The reference sequences and the mock sample are available at zenodo: https://doi.org/10.5281/zenodo.4588969

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
    │   ├── EMS_1.fq.gz
    │   ├── EMS_2.fq.gz
    │   ├── ENN_1.fq.gz
    │   ├── ENN_2.fq.gz
    │   ├── LHS_1.fq.gz
    │   ├── LHS_2.fq.gz
    │   ├── LNN_1.fq.gz
    │   └── LNN_2.fq.gz
    └── reference
        └── ZymoBIOMICS_genomes.fasta
        
This is already the expected input for LMAS. To execute LMAS with default resources and parameters you simply need to call the ``LMAS`` execution file either by typing:

.. code-block:: bash

    ./LMAS

Or alternatively:

.. code-block:: bash

    nextflow run main.nf
