# LMAS

![Nextflow CI](https://github.com/cimendes/LMAS/actions/workflows/ci_nextflow.yml/badge.svg)
![Pytest CI](https://github.com/cimendes/LMAS/actions/workflows/ci_templates.yml/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/lmas/badge/?version=latest)](https://lmas.readthedocs.io/en/latest/?badge=latest)

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.01.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with shifter](https://img.shields.io/badge/run%20with-shifter-lightgrey?labelColor=000000)](https://github.com/NERSC/shifter/)

[![DOI Dataset](https://zenodo.org/badge/DOI/10.5281/zenodo.4742651.svg)](https://doi.org/10.5281/zenodo.4742651)
[![Demo Report](https://img.shields.io/badge/Demo-lmas--demo.herokuapp.com%2F-blue)](https://lmas-demo.herokuapp.com/)


                      _    __  __   _   ___
       /\︵︵/\      | |  |  \/  | /_\ / __|
      (◕('人')◕)     | |__| |\/| |/ _ \\__ \
         |︶|        |____|_|  |_/_/ \_\___/

         Last Metagenomic Assembler Standing

## Overview

The *de novo* assembly of raw sequence data is a key process when analysing data from shotgun metagenomic sequencing. It allows recovering draft genomes from a pool of mixed raw reads, yielding longer sequences that offer contextual information and afford a more complete picture of the microbial community. It also represents one of the greatest bottlenecks when obtaining trustworthy, reproducible results.

**LMAS is an automated workflow enabling the benchmarking of traditional and metagenomic prokaryotic *de novo* assembly software using defined mock communities**. The results are presented in an interactive HTML report where selected global and reference specific performance metrics can be explored.

![LMAS Workflow](https://github.com/cimendes/LMAS/blob/main/docs/resources/LMAS_ECCMID.png)


## Instalation

Before installing LMAS, a few dependencies must be installed in your system:


### Nextflow

Nextflow (version 20.01.0 or higher) can be used on any POSIX compatible system (Linux, OS X, etc). It requires BASH and 
Java 8 (or higher) to be installed. More instructions are available [here](https://www.nextflow.io/docs/latest/getstarted.html).


### Container Engine

All components of LMAS are executed in docker containers, which means that you’ll need to have a container engine 
installed. The container engines available are the ones supported by Nextflow:

* [Docker](https://www.nextflow.io/docs/latest/docker.html),
* [Singularity](https://www.nextflow.io/docs/latest/singularity.html),
* [Shifter](https://github.com/NERSC/shifter)

If you already have any one of these installed, you are good to go as the provided docker containers are compatible 
with all engines available. If not, you’ll need to install one.


### Clone LMAS

You can clone this repository with `git clone git@github.com:cimendes/LMAS.git`, and all files will be in your local machine.


## Running LMAS

To use LMAS the following options are available:


                      _    __  __   _   ___
       /\︵︵/\      | |  |  \/  | /_\ / __|
      (◕('人')◕)     | |__| |\/| |/ _ \\__ \
         |︶|        |____|_|  |_/_/ \_\___/

         Last Metagenomic Assembler Standing

      Basic Usage: 
         nextflow run LMAS.nf

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

The reference sequences, in a single file, can be passed with the `--reference` parameter, and `--fastq` recieves the raw data for assembly.
The raw data is a collection of sequence fragments from the references, and can be either obtained *in silico* or from real
sequencing platforms. 

Users can customize the workflow execution either by using command line options or by modifying a simple plain-text 
configuration file (`configs/params.config`), where parameters are set as key-value pairs. The version of tools used can also 
be changed by providing new container tags in the appropriate configuration file (`contigs/containers.config`).

Users can select what profile to use with the `-profile` option. Several configurations are availabel in the profile configuration file (`contigs/profiles.config`). For a local execution we recommend running LMAS with either `-profile docker` or `-profile singularity`. HPC compatibility is available for SLURM, SGE, LSF, among others. 


## Output and Report

The output files are stored in the `results/` folder in the directory where the workflow was executed. 
The nextflow log file for the execution of the pipeline can be found in the directory of execution. Log files for each
of the components in the workflow are stored inside the `results/` folder.
LMAS creates an **interactive HTML report**, stored in the `report/` folder in the directory where the 
workflow was executed. To open the report simply click on the **index.html** file and the report will open on 
your default browser. 

LMAS comes pre-packaged with the JS source code for the interactive report, available in the `resources/` folder. 
The source code for the report is available in the [LMAS.js](https://github.com/B-UMMI/LMAS.js) repository. 

## Quick Start

### ZymoBIOMICS Microbial Community Standard

A bash script to download and structure the ZymoBIOMICS data to be used as input is provided (`get_data.sh`). 

      sh get_data.sh

Running this scipt downloads the [eight bacterial genomes and four plasmids of the ZymoBIOMICS Microbial Community Standards](https://zenodo.org/record/4588970#.YEeA83X7RhE) were used as reference. 
It contains tripled complete sequences for the following species:
- *Bacillus subtilis* 
- *Enterococcus faecalis*
- *Escherichia coli*
   - *Escherichia coli* plasmid
- *Lactobacillus fermentum*
- *Listeria monocytogenes*
- *Pseudomonas aeruginosa*
- *Salmonella enterica*
- *Staphylococcus aureus*
   - *Staphylococcus aureus* plasmid 1
   - *Staphylococcus aureus* plasmid 2
   - *Staphylococcus aureus* plasmid 3

It also downloads the raw sequence data of the mock communities, with an even ([ERR2984773](https://www.ebi.ac.uk/ena/browser/view/ERR2984773)) and logarithmic distribution of species ([ERR2935805](https://www.ebi.ac.uk/ena/browser/view/ERR2935805)), and the [complete reference sequences](https://zenodo.org/record/5579145/files/ZymoBIOMICS_genomes.fasta)

Simulated samples of the evenly and log distributed reads, with and without error, generated from the genomes in the Zymobiomics standard with [inSilicoSeq](https://github.com/HadrienG/InSilicoSeq) (version 1.5.2):
- ENN - Evenly distributed sample with no error model
- EMS - Evenly distributed sample with Illumina MiSeq error model
- LNN - Log distributed sample with no error model
- LHS - Log distributed sample with Illumina HiSeq error model

[DOI Dataset](https://doi.org/10.5281/zenodo.4588970)

After downloading the data you can simply run LMAS, with default parameters, with the following command:

      nextflow run LMAS.nf -profile docker

## Citation and Contacts

LMAS is developed at the Molecular [Microbiology and Infection Unit (UMMI)](http://darwin.phyloviz.net/wiki/doku.php) at the [Instituto de Medicina Molecular Joao Antunes](https://imm.medicina.ulisboa.pt/en/), in collaboration with [Microbiology, Advanced Genomics and Infection Control Applications Laboratory (MAGICAL)](https://morangiladlab.com) at the [Faculty of Health Sciences, Ben-Gurion University of the Negev](https://in.bgu.ac.il/en/fohs/Pages/default.aspx). 

This project is licensed under the [GPLv3 license](https://github.com/cimendes/LMAS/blob/main/LICENSE).

If you use LMAS please [cite this repository](https://github.com/cimendes/LMAS/blob/main/CITATION.cff).
