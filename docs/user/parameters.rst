Parameters
==========

A set of **default parameters** is provided, but these can be easily altered by either editing the 
``params.config`` file, or by passing the new value when executing the workflow with nextflow.
There are three main parameters in LMAS: **``--reference``, ``--fastq`` and ``--md``**. 

The **reference sequences**, in a single multifasta file, can be passed with the ``--reference`` parameter, and ``--fastq`` receives the 
**raw data** for assembly. The raw data is a collection of sequence fragments from the references and can be either 
obtained *in silico* or from real sequencing platforms. Users can pass text information, in a markdown file, 
on input samples to be presented in the report with the ``--md`` parameter.

Several options are available to alter the behaviour of the **assemblers** incorporated in LMAS, namely to alter 
the values of the k-mer for each assembly iteration. By default, these values reflect the corresponding default 
settings of the assemblers. 

The target values for some **quality assessment** metrics can also be adjusted, such as N50 and NG50.


Input Files
------------

fastq
^^^^^

Path expression to paired-end fastq files. Required.

* **Param:** :code:`--fastq`

* **Default:** :code:`data/fastq/*_{1,2}.*`


reference
^^^^^^^^^

Path to reference fasta file. Required.

* **Param:** :code:`--reference`

* **Default:** :code:`data/reference/*.fasta`


md
^^^

Path to markdown file with text to be displayed in the report. Optional.

* **Param:** :code:`--md`

* **Default:** :code:`data/*.md`


Assembly Quality Assessment
---------------------------

Minimum contig length
^^^^^^^^^^^^^^^^^^^^^
Value for minimum contig length, in basepairs.

* **Param:** :code:`--minLength`

* **Default:** 1000

Mapped reads threshold
^^^^^^^^^^^^^^^^^^^^^^^
Value for the minimum percentage of a read length aligning to the contig to be considered as mapped.

* **Param:** :code:`--mapped_reads_threshold`

* **Default:** 0.75

N Target
^^^^^^^^
Target value for the N*x*, NA*x* and NG*x* metrics. 

* **Param:** :code:`--n_target`

* **Default:** 0.9

L Target
^^^^^^^^
Target value for the L*x* metric. 

* **Param:** :code:`--l_target`

* **Default:** 0.5

L Target
^^^^^^^^
Scale of x-axis for the L, NA and NG metrics plots.

* **Param:** :code:`--plot_scale`

* **Default:** log

Assembler options
-----------------

ABySS
^^^^^^
* **Param:** :code:`--abyss`

* **Definition:** Boolean controling the execution of the ABySS assembler.

* **Default:** true

----------------

* **Param:** :code:`--abyssKmerSize`

* **Definition:** K-mer size for the ABySS assembler, as an intiger.

* **Default:** 96

-----------------

* **Param:** :code:`--abyssBloomSize`

* **Definition:** Bloom filter size for the ABySS assembler. It must be a sting with a value and an unit.

* **Default:** 2G

GATB Minia Pipeline
^^^^^^^^^^^^^^^^^^^

* **Param:** :code:`--gatb_minia`

* **Definition:** Boolean controling the execution of the GATB Minia Pipeline assembler.

* **Default:** true

----------------

* **Param:** :code:`--gatbKmerSize`

* **Definition:** K-mer sizes for the GATB Minia Pipeline assembler. It must be a sting with the values separated with a comma.

* **Default:** '21,61,101,141,181'

------------

* **Param:** :code:`--gatb_besst_iter`

* **Definition:** Number of iteration during BESST scaffolding for the GATB Minia Pipeline assembler.

* **Default:** 10000

------------

* **Param:** :code:`--gatb_error_correction`

* **Definition:** Boolean to control weather to skip error correction for the GATB Minia Pipeline assembler.

* **Default:** false

IDBA-UD
^^^^^^^^

* **Param:** :code:`--idba`

* **Definition:** Boolean controling the execution of the IDBA-UD assembler.

* **Default:** true

MetaHipMer2
^^^^^^^^^^^^^^^^^^^

* **Param:** :code:`--metahipmer2`

* **Definition:** Boolean controling the execution of the MetaHipMer2 assembler.

* **Default:** true

----------------

* **Param:** :code:`--metahipmer2KmerSize`

* **Definition:** K-mer sizes for the MetaHipMer2 assembler. It must be a sting with the values separated with a comma.

* **Default:** '21,33,55,77,99'


Minia
^^^^^

* **Param:** :code:`--minia`

* **Definition:** Boolean controling the execution of the minia assembler.

* **Default:** true

----------------

* **Param:** :code:`--miniaKmerSize`

* **Definition:** K-mer size for the minia assembler, as an intiger.

* **Default:** 31

MEGAHIT
^^^^^^^

* **Param:** :code:`--megahit`

* **Definition:** Boolean controling the execution of the MEGAHIT assembler.

* **Default:** true

----------------

* **Param:** ``--megahitKmerSize``

* **Definition:** K-mer sizes for the MEGAHIT assembler. It must be a sting with the values separated with a comma.

* **Default:** '21,29,39,59,79,99,119,141'

metaSPAdes
^^^^^^^^^^

* **Param:** :code:`--metaspades`

* **Definition:** Boolean controling the execution of the metaSPAdes assembler.

* **Default:** true

----------------

* **Param:** :code:`--metaspadesKmerSize`

* **Definition:** K-mer sizes for the metaSPAdes assembler. It must be a sting with 'auto' or with the values separated with a space.

* **Default:** 'auto'

SPAdes
^^^^^^
* **Param:** :code:`--spades`

* **Definition:** Boolean controling the execution of the SPAdes assembler.

* **Default:** true

----------------

* **Param:** :code:`--spadesKmerSize`

* **Definition:** K-mer sizes for the metaSPAdes assembler. It must be a sting with 'auto' or with the values separated with a space.

* **Default:** 'auto'

SKESA
^^^^^^^^

* **Param:** :code:`--skesa`

* **Definition:** Boolean controling the execution of the SKESA assembler.

* **Default:** true

Unicycler
^^^^^^^^^^

* **Param:** :code:`--unicycler`

* **Definition:** Boolean controling the execution of the Unicycler assembler.

* **Default:** true

VelvetOptimiser
^^^^^^^^^^^^^^^

* **Param:** :code:`--velvetoptimiser`

* **Definition:** Boolean controling the execution of the VelvetOptimiser assembler.

* **Default:** true

------------  

* **Param:** :code:`--velvetoptimiser_hashs`

* **Definition:** Starting K-mer size for the VelvetOptimiser assembler, as an intiger.

* **Default:** 19

------------  

* **Param:** :code:`--velvetoptimiser_hashe`

* **Definition:** End K-mer size for the VelvetOptimiser assembler, as an intiger.

* **Default:** 31


Execution Resources Parameters
-------------------------------

CPUs
^^^^^^^^
Number of CPUs for the assembly and mapping processes, as an intiger.
This resource is double for each retry until max_cpus is reached.

* **Param:** :code:`--cpus`

* **Default:** 8

Memory
^^^^^^^^
Memory for the assembly and mapping processes, in the format of 'value'.'unit'.
This resource is double for each retry until max_memory is reached.

* **Param:** :code:`--memory`

* **Default:** 32.Gb

Time
^^^^^^^^
Time limit for the assembly and mapping processes, in the format of 'value'.'unit'.
This resource is double for each retry until max_time is reached.

* **Param:** :code:`--memory`

* **Default:** 24.h 

Max_cpus
^^^^^^^^
Maximum number of CPUs for the assembly and mapping processes, as an intiger.

* **Param:** :code:`--max_cpus`

* **Default:** 32

Max_memory
^^^^^^^^^^^^
Maximum memory for the assembly and mapping processes, in the format of 'value'.'unit'.

* **Param:** :code:`--max_memory`

* **Default:** 100.Gb

Max_time
^^^^^^^^^^^^
Maximum time for the assembly and mapping processes, in the format of 'value'.'unit'.

* **Param:** :code:`--max_memory`

* **Default:** 72.h