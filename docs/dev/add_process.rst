Add Assembler Process
=====================

New assemblers can be added with minimal changes to the pipeline, 
so that LMAS can be expanded as novel algorithms are developed. It's implementation in DSL2 
greatly facilitates this process. 

The assemblers implemented are available in the **assembly** module located in the  
`modules folder <https://github.com/B-UMMI/LMAS/tree/main/modules/assembly>`_.
The `assembly.nf <https://github.com/B-UMMI/LMAS/tree/main/modules/assembly/assembly.nf>`_ is the nextflow file that contains all the assembly processes.

The current available assemblers are:

* `ABySS <https://github.com/bcgsc/abyss>`_
* `BCALM 2 <https://github.com/GATB/bcalm>`_ 
* `GATB-Minia Pipeline <https://github.com/GATB/gatb-minia-pipeline>`_
* `IDBA <https://github.com/loneknightpy/idba>`_
* `MEGAHIT <https://github.com/voutcn/megahit>`_
* `MetaHipMer2 <https://bitbucket.org/berkeleylab/mhm2>`_
* `metaSPAdes <https://github.com/ablab/spades>`_
* `MINIA <https://github.com/GATB/minia>`_
* `SKESA <https://github.com/ncbi/SKESA>`_
* `SPAdes <https://github.com/ablab/spades>`_
* `Unicycler <https://github.com/rrwick/Unicycler>`_
* `VelvetOptimiser <https://github.com/tseemann/VelvetOptimiser>`_

Detailed information is available in the `Short-Read (Meta)Genomic Assemblers <../user/assemblers.html>`_ page.

.. warning:: To add an assembler, it must be ensured that **short-read paired-end sequence data** can be provided as input. 


Changing assembler version
-----------------------------------

The easiest way to change a version of a particular assembler in LMAS is by changing the containers for the assembler process.
This is done through altering the container property in the `containers.config <https://github.com/B-UMMI/LMAS/blob/main/conf/containers.config>`_ file.

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
    in the container for it to work with LMAS or any other Nextflow workflow.

Adding a new assembler
-----------------------------------

Create an issue with an assembler suggestion
:::::::::::::::::::::::::::::::::::::::::::::::

An issue template is available to collect the necessary information for an assembler to be added to LMAS.
Some information is required:

* Container for the execution of the assembler, containing the executable in the PATH and Nextflow's ps dependency;
* Command to capture the assembler version, if available;
* Minimal command to execute the assembler with short-read paired-end sequencing datasets;
* Parameters (such as k-mer lists) to be passed onto the assembler.

By default, all assemblies are run with 8 CPUs and 32GB of memory. 


Add process to assembly.nf manually
:::::::::::::::::::::::::::::::::::::::::::::

To add a new assembler to LMAS, a few steps must be completed. All alterations 
needed will be perfomed in the `assembly.nf <https://github.com/B-UMMI/LMAS/tree/main/modules/assembly/assembly.nf>`_ file,
the `params.config <https://github.com/B-UMMI/LMAS/tree/main/conf/params.config>`_ file and the
`containers.config <https://github.com/B-UMMI/LMAS/tree/main/conf/containers.config>`_.

1. Add the needed parameters

In the the `params.config <https://github.com/B-UMMI/LMAS/tree/main/conf/params.config>`_ file,
add a new key-value pair for any parameter necessary to run the assembler, such as the list of k-mer values to use. 
The fastq input data is passed through the main `--fastq` parameter so it should not be included.

.. warning:: All assemblers in LMAS are toggleble through a `--<assembler name>` parameter, and this should be included in this file.

2. **Add a new process with the assembler**

In the `assembly.nf <https://github.com/B-UMMI/LMAS/tree/main/modules/assembly/assembly.nf>`_ file, you need
to add the process to execute the new assembler in the section marked with `\\PROCESSES`. 

To create the new process, you can use the following template, substituting ``NEW_ASSEMBLER`` with the new
assembler name:

.. code-block:: bash

    process NEW_ASSEMBLER {
        tag { sample_id }
        label 'process_assembly'
        publishDir 'results/assembly/NEW_ASSEMBLER/'

        input:
        set sample_id, path(fastq) 
        val kmers from IN_NEW_ASSEMBLER_kmers

        when:
        pararm.NEW_ASSEMBLER

        output:
        set sample_id, val("NEW_ASSEMBLER"), file('*.fasta'), emit: assembly
        file(".*version"), emit: version

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

You can access this values in the ``.nf`` file with ``params.<parameter>``.
For example:

.. code-block:: bash

    IN_NEW_ASSEMBLER_kmers = Channel.value(params.newassemblerKmers)

.. warning:: Parameters need to be passed into a process through a channel. 

This should be added inside the `assembly_wf` worflow in the end of the file.

Additionally, The new process needs to be added in the `main:` section of the
workflow. 

3. **Add assembly to main assembly collection**

The channel with the version information must be merged into the main assembly collection channel, emitted by the `assembly_wf` workflow.

It should look like:

.. code-block:: bash

    all_assemblies = ABYSS.out.assembly | mix(BCALM2.out.assembly, 
                                              GATBMINIAPIPELINE.out.assembly,
                                              IDBA.out.assembly,
                                              MEGAHIT.out.assembly,
                                              METAHIPMER2.out.assembly,
                                              METASPADES.out.assembly,
                                              MINIA.out.assembly,
                                              NEW_ASSEMBLER.out.version, // new channel added 
                                              SKESA.out.assembly,
                                              SPADES.out.assembly,
                                              UNICYCLER.out.assembly,
                                              VELVETOPTIMISER.out.assembly)

.. warning:: To facilitate reading, please respect the alphabetical order. 

4. **Add version to main version collection**

The channel with the version information must be merged into the main version collection channel, emitted by the `assembly_wf` workflow.

It should look like:

.. code-block:: bash

    all_versions = ABYSS.out.version | mix(BCALM2.out.version, 
                                           GATBMINIAPIPELINE.out.version,
                                           IDBA.out.version,
                                           MEGAHIT.out.version,
                                           METAHIPMER2.out.version,
                                           METASPADES.out.version,
                                           MINIA.out.version,
                                           NEW_ASSEMBLER.out.version,  // new channel added 
                                           SKESA.out.version,
                                           SPADES.out.version,
                                           UNICYCLER.out.version,
                                           VELVETOPTIMISER.out.version) | collect

.. warning:: To facilitate reading, please respect the alphabetical order. 


5. **Add the container for the new assembler**

The container for the new assembler need to be added to the ``container.config`` file
in the `conf/` directory.

It should look like:

.. code-block:: bash
    
    withName: NEW_ASSEMBLER {
        container = "<repository>/NEW_ASSEMBLER:<tag>"
    }
