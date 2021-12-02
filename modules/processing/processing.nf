nextflow.enable.dsl=2

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
    template "process_reads.py"
}

process PROCESS_VERSION {

    label 'process_script'

    input:
    file version

    output:
    file 'versions.json'

    script:
    template "process_versions.py"
}