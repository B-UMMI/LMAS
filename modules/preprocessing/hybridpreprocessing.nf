nextflow.enable.dsl=2

// PROCESSES
process PROCESS_REFERENCE {

    label 'process_script'

    input:
    path reference_fasta 

    output:
    path('triple_reference.fasta')

    script:
    template "process_reference.py"
}

process PROCESS_READS {

    tag { sample_id }
    label 'process_script'

    input:
    tuple val(sample_id), path(fastq) 

    output:
    path('*_reads_report.json')

    script:
    template "process_reads_hybrid.py"
}

// WORKFLOW
workflow hybridpreprocessing_wf {

    take:
    reference
    fastq_short
    fastq_log

    main:
    PROCESS_REFERENCE(reference)
    PROCESS_READS(fastq_short, fastq_log)

    emit:
    triple_reference = PROCESS_REFERENCE.out
    reads_info = PROCESS_READS.out
}