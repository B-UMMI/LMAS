Add Assembler Process
=====================

New assemblers can be added with minimal changes to the pipeline, 
so that LMAS can be expanded as novel algorithms are developed.

The assemblers implemented are available in the **main file** in 
`LMAS.nf <https://github.com/cimendes/LMAS/blob/main/LMAS.nf>`_, to be executed by Nextflow.

The current available assemblers are:
* BCALM2
* GATB-Minia Pipeline
* IDBA
* MINIA
* MEGAHIT
* METASPADES
* SKESA
* SPADES
* UNICYCLER
* VELVETOPTIMIZER

To add an assembler, it must be insured that **short-read paired-end sequence data** can be provided as input. 

More information is available at `assemblers <../user/assemblers>`_.

Changing assembler version
-----------------------------------

The easiest way to change a version of a particular assembler in LMAS is by changing the containers for the assembler process.
This is done though altering the container property in the `containers.config <https://github.com/cimendes/LMAS/blob/main/containers.config>`_ file.

For example, for the ``SPADES`` process, the container "cimendes/spades:3.15.0-1" can be altered to another one that implements a
different version of the tool: 

.. code-block:: bash
    withName: SPADES {
                container = "cimendes/spades:3.15.0-1"
            }


.. code-block:: bash
    withName: SPADES {
                container = "cimendes/spades:3.14.1-1"
            }

.. warning:: You must ensure that the assembler executable is available in the $PATH and that ps is installed 
    in the container for it to work with LMAS.

Adding a new assembler
-----------------------------------

Create an issue with an assembler suggestion
:::::::::::::::::::::::::::::::::::::::::::::::

An issue template is available to collect the necessary information for an assembler to be added to LMAS.
Some information is required:

* Container for the execution of the assembler, containing the executable in the PATH and Nextflow's ps dependency
* Command to capture the assembler version, if available,
* Minimal command to execute the assembler with short-read paired-end sequencing datasets
* Parameters (such as k-mer lists) to be passed onto the assembler

By default, all assemblies are run with 4 CPUS and 16Gb of memory. 


Add process to LMAS.nf manually
:::::::::::::::::::::::::::::::::
To add a new assembler to LMAS, a few steps must be completed:

1. **Add new channel for the FASTQ datasets**

A new channel to provide the raw sequence data to the new assembler must be created.
Simply add a new channel, named for example ``IN_<NEW_ASSEMBLER>``, to the ``into`` operator
that splits the data in the ``IN_fastq_raw`` channel in this `line <https://github.com/cimendes/LMAS/blob/main/LMAS.nf#L58>`_.

It should look like:

.. code-block:: bash
    // SET CHANNELS FOR ASSEMBLERS
    IN_fastq_raw.into{
        IN_PROCESS_READS;
        IN_BCALM2;
        IN_GATB_MINIA_PIPELINE;
        IN_MINIA;
        IN_MEGAHIT;
        IN_METASPADES;
        IN_UNICYCLER;
        IN_IDBA;
        IN_SPADES;
        IN_SKESA;
        IN_VELVETOPTIMIZER;
        IN_NEW_ASSEMBLER; // new channel added
        IN_TO_MAP} //mapping channel - minimap2

.. warning:: Make sure the channel name isn't used elsewere. Otherwise Nextflow will throw an error

2. **Add new process with the assembler**

Parameters to be passed on to this new process can be added in the `params.config <https://github.com/cimendes/LMAS/blob/main/params.config>`_ file.
You can access this values in the ``.nf`` file with ``params.<parameter>``.
For example

.. code-block:: bash
    IN_NEW_ASSEMBLER_kmers = Channel.value(params.newassemblerKmers)

.. warning:: For parameters need to be passed into a process through a channel. 

To create the new process, you can use the following template, substituting ``NEW_ASSEMBLE`` with the new
assembler name:

.. code-block:: bash

    process NEW_ASSEMBLER {
        tag { sample_id }
        publishDir 'results/assembly/NEW_ASSEMBLER/'

        input:
        set sample_id, file(fastq_pair) from IN_NEW_ASSEMBLER
        val kmers from IN_NEW_ASSEMBLER_kmers

        output:
        set sample_id, val("NEW_ASSEMBLER"), file('*.fasta') into OUT_NEW_ASSEMBLER
        file(".*version") into NEW_ASSEMBLER_VERSION

        script:
        """
        // capture assembler version and save into 
        <version command> > .${sample_id}_NEWASSEMBLER_version

        // Run assembly in a try-except 
        {
            <assembly command>
            echo pass > .status
        } || {
            echo fail > .status
        }
        """
    }

.. warning:: You can access each of the fastq files with ${fastq_pair[1]} and ${fastq_pair[2]}.


3. **Add version to main version collection**

The channel with the version information must be merged into the main version collection channel
for it to be processed accordingly in this `line <https://github.com/cimendes/LMAS/blob/main/LMAS.nf#L422>`_.

It should look like:

.. code-block:: bash
    // VERSION COLLECTION
    BCALM2_VERSION.mix(GATB_VERSION,
                        MINIA_VERSION,
                        MEGAHIT_VERSION,
                        METASPADES_VERSION,
                        UNICYCLER_VERSION,
                        SPADES_VERSION,
                        SKESA_VERSION,
                        VELVETOPTIMIZER_VERSION,
                        NEW_ASSEMBLER_VERSION,  // new channel added 
                        IDBA_VERSION).set{ALL_VERSIONS}

4. **Add assembly to main assembly collection**

The channel with the assembly produced  must be merged into the main assembly collection channel
for it to be processed. This is done in this `line <https://github.com/cimendes/LMAS/blob/main/LMAS.nf#L445>`_.

It should look like:

.. code-block:: bash
    // ASSEMBLY COLLECTION
    OUT_BCALM2.mix(OUT_GATB,
                    OUT_MINIA,
                    OUT_MEGAHIT,
                    OUT_METASPADES,
                    OUT_UNICYCLER,
                    OUT_SPADES,
                    OUT_SKESA,
                    OUT_VELVETOPTIMIZER,
                    OUT_NEW_ASSEMBLER,   // new channel added 
                    OUT_IDBA).set{ALL_ASSEMBLERS}