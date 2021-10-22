process PROCESS_REFERENCE {

    label 'process_script'

    input:
    file reference_fasta from TO_TRIPLE

    output:
    file 'triple_reference.fasta' into OUT_REFERENCE_TRIPLE

    script:
    template "process_reference.py"
}

process PROCESS_READS {

    tag { sample_id }
    label 'process_script'

    input:
    tuple sample_id, file(fastq) from IN_PROCESS_READS

    output:
    file '*_reads_report.json' into PROCESS_READS

    script:
    template "process_reads.py"
}