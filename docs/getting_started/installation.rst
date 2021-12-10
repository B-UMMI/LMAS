Installation
============

LMAS can be installed through Github (https://github.com/cimendes/LMAS).
It requires a `Nextflow <https://www.nextflow.io/>`_ installation (version ≥ 21.04.1) 
and can be used on any POSIX compatible system (Linux, OS X, etc). All components of LMAS are executed in `Docker containers <https://www.docker.com/>`_, 
being a container engine required. 

Nextflow allows integration with multiple alternatives, such as `Shifter <https://github.com/NERSC/shifter/>`_ or 
`Singularity <https://singularity.hpcng.org/>`_, so a particular one isn’t required. 

To ensure the robustness of LMAS workflow and the custom python code for the quality assessment of assemblies, **continuous integration** of both the main workflow
and the python templates is performed with `GitHub Actions <https://github.com/features/actions>`_ and `pytest <https://docs.pytest.org/en/6.2.x/>`_. 

Below it's a step by step guide on how to install LMAS and all its dependencies.

Step 1. Nextflow
-----------------

`Nextflow <https://www.nextflow.io/>`_ (version 20.01.0 or higher) can be used on any POSIX compatible system (Linux, OS X, etc). 
It requires BASH and Java 8 (or higher) to be installed. 

.. important::

    Instructions on how to install Nextflow are available `here <https://www.nextflow.io/docs/latest/getstarted.html>`_

Step 2. Container engine
-------------------------

All components of LMAS are executed in docker containers, which means that you’ll need to have a container engine 
installed. The container engines available are the ones supported by Nextflow:

- `Docker`_,
- `Singularity`_
- Shifter (undocumented)

If you already have any one of these installed, you are good to go as the provided docker containers are compatible 
with all engines available. If not, you’ll need to install one.


Singularity
:::::::::::

Singularity is available to download and install `here <http://singularity.lbl.gov/install-linux>`_.
Make sure that you have singularity v2.5.x or higher.
Note that singularity should be installed as root and available on the machine(s) that
will be running the nextflow processes.

.. important::

    Singularity is available as a bioconda package. However, conda installs singularity
    in user space without root privileges, which may prevent singularity images from
    being correctly downloaded. **Therefore it is not recommended that you install
    singularity via bioconda**.

Docker
::::::

Docker can be installed following the instructions on the website:
https://www.docker.com/community-edition#/download.
To run docker as a non-root user, you'll need to follow the instructions
on the website: https://docs.docker.com/install/linux/linux-postinstall/#manage-docker-as-a-non-root-user

Step 3. Clone LMAS
-------------------

You can clone this repository with git.

.. code-block:: bash

    git clone https://github.com/cimendes/LMAS.git 

All files will be on your local machine.

To run LMAS you can simply call it with:

.. code-block:: bash

    ./LMAS <options>
   
If no option or `--help` is provided, LMAS will display its help message. Otherwise, the `--fastq` and `--reference` options are mandatory. By default they are set to `'data/fastq/*_{1,2}.*'` and `'data/reference/*.fasta'` respectively.

The main execution file for Nextflow is ``main.nf``. Alternatively you can call LMAS directly with Nextflow:

.. code-block:: bash

      nextflow run main.nf <options>

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

