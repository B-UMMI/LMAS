nextflow.enable.dsl=2

// PROCESSES
process METASPADES {

    tag { sample_id }
    label 'process_assembly'
    publishDir "results/$sample_id/assembly/metaSPAdes/"

    when:
    params.metaspades

    input:
    tuple val(sample_id), path(fastq_pair)
    tuple val(ont_id), path(fastq_ont)
    val kmers

    output:
    tuple val(sample_id), val('metaSPAdes'), path('*_metaspades.fasta'), emit: assembly
    path('.*version'), emit: version

    script:
    """
    metaspades.py --version &> version
    cat version | awk -F ' ' '{print \$4}' | awk -F 'v' '{print \$2}' > .${sample_id}_metaSPAdes_version
    rm version
    {
        metaspades.py --only-assembler --threads $task.cpus -k $kmers \
        -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} --nanopore ${fastq_ont}\
        -o metaspades
 
        mv metaspades/contigs.fasta ${sample_id}_metaspades.fasta
        echo pass > .status
    } || {
        echo fail > .status
        :> ${sample_id}_metaspades.fasta
    }
    rm -r metaspades || true
    """
}

process SPADES {

    tag { sample_id }
    label 'process_assembly'
    publishDir "results/$sample_id/assembly/SPAdes/", pattern: '*.fasta'

    when:
    params.spades

    input:
    tuple val(sample_id), file(fastq_pair) 
    tuple val(ont_id), path(fastq_ont)
    val kmers

    output:
    tuple val(sample_id), val('SPAdes'), path('*_spades.fasta'), emit: assembly
    path('.*version'), emit: version

    script:
    """
    spades.py --version &> version
    cat version | awk -F ' ' '{print \$4}' | awk -F 'v' '{print \$2}' > .${sample_id}_SPAdes_version
    rm version
    {
        spades.py --only-assembler --threads $task.cpus -k $kmers \
        -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} --nanopore ${fastq_ont}\
        -o spades

        mv spades/contigs.fasta ${sample_id}_spades.fasta
    } || {
        echo fail > .status
        :> ${sample_id}_spades.fasta
    }
    rm -r spades || true
    """
}

process UNICYCLER {

    tag { sample_id }
    label 'process_assembly'
    publishDir "results/$sample_id/assembly/unicycler"

    when:
    params.unicycler

    input:
    tuple val(sample_id), file(fastq_pair)
    tuple val(ont_id), path(fastq_ont)

    output:
    tuple val(sample_id), val('Unicycler'), path('*_unicycler.fasta'), emit: assembly
    path('.*version'), emit: version

    script:
    """
    unicycler --version | awk -F ' v' '{print \$2}' | awk NF > .${sample_id}_Unicycler_version
    {
        unicycler -t $task.cpus -o . --no_correct --no_pilon \
        -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} -l ${fastq_ont}

        mv assembly.fasta ${sample_id}_unicycler.fasta
        echo pass > .status
    } || {
        echo fail > .status
        :> ${sample_id}_unicycler.fasta
    }
    rm *best_spades_graph* *overlaps_removed* *bridges_applied* *final_clean* || true
    """
}

// WORKFLOWS
workflow hybridassembly_wf {

    metaspadesKmerSize = Channel.value(params.metaspadesKmerSize)
    spadesKmerSize = Channel.value(params.spadesKmerSize)

    take:
    IN_fastq_raw
    IN_fastq_ont

    main:
    METASPADES(IN_fastq_raw, IN_fastq_ont, metaspadesKmerSize)
    SPADES(IN_fastq_raw, IN_fastq_ont, spadesKmerSize)
    UNICYCLER(IN_fastq_raw, IN_fastq_ont)


    emit:
    all_assemblies = METASPADES.out.assembly | mix(SPADES.out.assembly,
                                                   UNICYCLER.out.assembly)
    all_versions = METASPADES.out.version | mix(SPADES.out.version,
                                                UNICYCLER.out.version) | collect
}