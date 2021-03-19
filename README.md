# LMAS

[![Documentation Status](https://readthedocs.org/projects/lmas/badge/?version=latest)](https://lmas.readthedocs.io/en/latest/?badge=latest)
[![DOI Dataset](https://zenodo.org/badge/DOI/10.5281/zenodo.4588970.svg)](https://doi.org/10.5281/zenodo.4588970)



                      _    __  __   _   ___
       /\︵︵/\      | |  |  \/  | /_\ / __|
      (◕('人')◕)     | |__| |\/| |/ _ \\__ \
         |︶|        |____|_|  |_/_/ \_\___/

         Last Metagenomic Assembler Standing

## Table of Contents

  * [Overview](#overview)
  * [Instalation](#instalation)
    + [Nextflow](#nextflow)
    + [Container Engine](#container-engine)
    + [Clone LMAS](#clone-lmas)
  * [Running LMAS](#running-lmas)
  * [Customizing LMAS](#customizing-lmas)
  * [Output and Report](#output-and-report)
  * [Proof of concept](#proof-of-concept)
    + [ZymoBIOMICS Microbial Community Standard](#zymobiomics-microbial-community-standard)
  * [Citation and Contacts](#citation-and-contacts)


## Overview

The *de novo* assembly of raw sequence data is a key process when analysing data from shotgun metagenomic sequencing. It allows recovering draft genomes from a pool of mixed raw reads, yielding longer sequences that offer contextual information and afford a more complete picture of the microbial community. It also represents one of the greatest bottlenecks when obtaining trustworthy, reproducible results.

LMAS is an automated workflow enabling the benchmarking of traditional and metagenomic prokaryotic *de novo* assembly software using defined mock communities. The results are presented in an interactive HTML report where selected global and reference specific performance metrics can be explored.

LMAS expects **tripled reference sequences**. Each reference should be provided in a single FASTA file where the linearized 
reference replicons are concatenated three times to ensure that contigs can fully align even with start-end overlap and 
regardless of their starting position relative to that of the reference. Read data, in **paired-end form*, reflecting the 
sequences in the reference genomes, is required to be passed on for assembly and downstream quality assessment. 

The mock communities can be provided by the user to better reflect the samples of interest. New assemblers can be added with minimal changes to the pipeline, so that LMAS can be expanded as novel algorithms are developed.

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
* [Shifter](https://github.com/NERSC/shifter) (undocumented)

If you already have any one of these installed, you are good to go as the provided docker containers are compatible 
with all engines available. If not, you’ll need to install one.

### Clone LMAS

You can clone this repository with `git clone https://github.com/cimendes/LMAS.git`, and all files will be in your local machine.

## Running LMAS

After you have a local installation of LMAS, the mock community data needs to be downloaded.

The triple reference sequences can be passed with the `--reference` parameter, and `--fastq` recieves the raw data for assembly.
The raw data is a collection of sequence fragments from the references, and can be either obtained *in silico* or from real
sequencing platforms. 


## Customizing LMAS

Users can customize the workflow execution either by using command line options or by modifying a simple plain-text 
configuration file (`params.config`), where parameters are set as key-value pairs. The version of tools used can also 
be changed by providing new container tags in the appropriate configuration file (`containers.config`), as well as the 
resources for each process (`resources.config`).


## Output and Report

The output files are stored in the `results/` folder in the directory where the workflow was executed. 
The nextflow log file for the execution of the pipeline can be found in the directory of execution. Log files for each
of the components in the workflow are stored inside the `results/` folder.
LMAS creates an **interactive HTML report**, stored in the `report/` folder in the directory where the 
workflow was executed. To open the report simply click on the **index.html** file and the report will open on 
your default browser. 

LMAS comes pre-packaged with the JS source code for the interactive report, available in the `resources/` folder. 
The source code for the report is available in the [lmas_report](https://github.com/cimendes/lmas_report) repository. 

## Evaluation

### ZymoBIOMICS Microbial Community Standard

A script to download and structure the ZymoBIOMICS data to be used as input is provided (`get_data.sh`). 

Running this scipt downloads the [eight bacterial genomes and four plasmids of the ZymoBIOMICS Microbial Community Standards](https://zenodo.org/record/4588970#.YEeA83X7RhE) were used as the triple reference. 
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

It also downloads the raw sequence data of the mock communities, with an even and logarithmic distribution of species ([ERR2984773](https://www.ebi.ac.uk/ena/browser/view/ERR2984773) and [ERR2935805](https://www.ebi.ac.uk/ena/browser/view/ERR2935805)), and a simulated sample of the evenly distributed reads generated from the genomes in the Zymobiomics standard ([mockSample](https://zenodo.org/record/4588970#.YEeA83X7RhE)).

The resulting LMAS report is available at [https://lmas-demo.herokuapp.com](https://lmas-demo.herokuapp.com)

## Citation and Contacts

LMAS is developed at the Molecular [Microbiology and Infection Unit (UMMI)](http://darwin.phyloviz.net/wiki/doku.php) at the [Instituto de Medicina Molecular Joao Antunes](https://imm.medicina.ulisboa.pt/en/), in collaboration with [Microbiology, Advanced Genomics and Infection Control Applications Laboratory (MAGICAL)](https://morangiladlab.com) at the [Faculty of Health Sciences, Ben-Gurion University of the Negev](https://in.bgu.ac.il/en/fohs/Pages/default.aspx). 

This project is licensed under the [GPLv3 license](https://github.com/cimendes/LMAS/blob/main/LICENSE).

If you use LMAS please cite this repository.
