Create a Mock Community for LMAS
==================================

Triple Reference
----------------

.. warning:: In construction

Read Data
---------

For the reaction of simulated sequencing data, we recommend using `InSilicoSeq <https://github.com/HadrienG/InSilicoSeq>`_, 
a simulator that produces realistic Illumina reads primarily intended for simulating metagenomic samples, but that can also 
be used to produce sequencing data from a single genome.

It provides models for Illumina HiSeq, NovaSeq and Miseq to realistically estimate the read quality of real sequencing data.

Installation
:::::::::::::

`InSilicoSeq <https://github.com/HadrienG/InSilicoSeq>`_ can be installed througn conda or pip. It requires ``python >= 3.5``. 
.. code-block:: bash
conda install -c bioconda insilicoseq
pip install insilicoseq

Alternatively, a `docker container <https://hub.docker.com/r/hadrieng/insilicoseq>`_ is available.

Usage
::::::

To generate mock communities, the following command can be used:

.. code-block:: bash
iss generate --genomes triplegenomes.fasta --model hiseq --output issreads 

where ``triplegenomes.fasta`` should be replaced by a (multi-)fasta file containing the triple reference genome(s) 
from which the simulated reads will be generated.

InSilicoSeq comes with 3 error models that can be passed with the ``--model`` parameter: ``MiSeq``, ``HiSeq`` and ``NovaSeq``
You can change the number of cpus with the ``--cpus`` parameters, and the number of reads to generate with the 
``-n`` parameter. 

The mock reads will be saved with the ``issreads`` prefix. 




