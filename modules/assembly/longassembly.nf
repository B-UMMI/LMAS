nextflow.enable.dsl=2

// PROCESSES
process CANU {

    tag { sample_id }
    label 'process_long_assembly'
    publishDir "results/$sample_id/assembly/CANU/"

    when:
    params.canu

    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), val('CANU'), path('*_CANU.fasta'), emit: assembly
    path('.*version'), emit: version

    script:
    """
    canu --version | awk -F ' ' '{print \$2}' > .${sample_id}_CANU_version
    {
        canu -p assembly -d canu_out \
        genomeSize="${params.canu_genomesize}" -nanopore "${fastq}" \
        useGrid="false" minInputCoverage=0.1 stopOnLowCoverage=0.1
        
        mv canu_out/assembly.contigs.fasta ${sample_id}_CANU.fasta
        echo pass > .status
    } || {
        echo fail > .status
        :> ${sample_id}_CANU.fasta
    }
    rm -r canu_out || true
    """
}

process RAVEN {

    tag { sample_id }
    label 'process_long_assembly'
    publishDir "results/$sample_id/assembly/RAVEN/"

    when:
    params.raven

    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), val('RAVEN'), path('*_RAVEN.fasta'), emit: assembly
    path('.*version'), emit: version

    script:
    """
    raven --version > .${sample_id}_RAVEN_version
    {
        raven -t $task.cpus ${fastq} > ${sample_id}_RAVEN.fasta
        echo pass > .status
    } || {
        echo fail > .status
        :> ${sample_id}_RAVEN.fasta
    }
    rm -r raven.cereal || true
    """
}

process FLYE {

    tag { sample_id }
    label 'process_long_assembly'
    publishDir "results/$sample_id/assembly/FLYE/"

    when:
    params.flye

    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), val("FLYE"), path("*_FLYE.fasta"), emit: assembly
    path(".*version"), emit: version

    script:
    """
    flye --version > .${sample_id}_FLYE_version
    {
        flye --nano-raw ${fastq} --out-dir flye_out -t $task.cpus
        mv flye_out/assembly.fasta ${sample_id}_FLYE.fasta
        echo pass > .status
    } || {
        echo fail > .status
        :> ${sample_id}_FLYE.fasta
    }
    rm -r flye_out || true
    """
}

process METAFLYE {

    tag {sample_id}
    label 'process_long_assembly'
    publishDir "results/$sample_id/assembly/METAFLYE/"

    when:
    params.metaflye

    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), val("METAFLYE"), path("*_METAFLYE.fasta"), emit: assembly
    path(".*version"), emit: version

    script:
    """
    flye --version > .${sample_id}_METAFLYE_version
    {
        flye --nano-raw ${fastq} --meta --out-dir flye_out -t $task.cpus
        mv flye_out/assembly.fasta ${sample_id}_METAFLYE.fasta
        echo pass > .status
    } || {
        echo fail > .status
        :> ${sample_id}_METAFLYE.fasta
    }
    rm -r flye_out || true
    """
}

process RA{

    tag {sample_id}
    label 'process_long_assembly'
    publishDir "results/$sample_id/assembly/ra/"

    when:
    params.ra

    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), val("RA"), path("*_RA.fasta"), emit: assembly
    path(".*version"), emit: version

    script:
    """
    ra --version | awk -F 'v' '{print \$2}' | awk NF > .${sample_id}_RA_version
    {
        ra -t $task.cpus -x ont ${fastq[0]} > ${sample_id}_RA.fasta
        echo pass > .status
    } || {
        echo fail > .status
        :> ${sample_id}_RA.fasta
    }
    """
}

process WTDBG2{

    tag {sample_id}
    label 'process_long_assembly'
    publishDir "results/$sample_id/assembly/wtdbg2/"

    when:
    params.wtdbg2

    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), val("WTDBG2"), path("*_WTDBG2.fasta"), emit: assembly
    path(".*version"), emit: version

    script:
    """
    wtdbg2 -V | awk -F ' ' '{print \$2}' | awk NF > .${sample_id}_WTDBG2_version
    {
        wtdbg2 -t $task.cpus -i ${fastq[0]} -o ${sample_id}_WTDBG2
        wtpoa-cns -t $task.cpus -i ${sample_id}_WTDBG2.ctg.lay.gz -fo ${sample_id}_WTDBG2.fasta
        echo pass > .status
    } || {
        echo fail > .status
        :> ${sample_id}_WTDBG2.fasta
    }
    """
}


// WORKFLOWS
workflow longassembly_wf {

    take:
    IN_fastq_raw

    main:
    CANU(IN_fastq_raw)
    RAVEN(IN_fastq_raw)
    FLYE(IN_fastq_raw)
    METAFLYE(IN_fastq_raw)
    RA(IN_fastq_raw)
    WTDBG2(IN_fastq_raw)

    emit:
    all_assemblies = CANU.out.assembly | mix(RAVEN.out.assembly,
                                              FLYE.out.assembly,
                                              METAFLYE.out.assembly,
                                              RA.out.assembly,
                                              WTDBG2.out.assembly)
    all_versions = CANU.out.version | mix(RAVEN.out.version,
                                           FLYE.out.version,
                                           METAFLYE.out.version,
                                           RA.out.version,
                                           WTDBG2.out.version) | collect
}
