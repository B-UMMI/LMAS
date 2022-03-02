/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { preprocessing_wf } from '../modules/preprocessing/preprocessing'
include { assembly_wf } from '../modules/assembly/longassembly'

/*
========================================================================================
    RUN WORKFLOW
========================================================================================
*/

workflow LONGLMAS {

    take:
    IN_reference_raw
    IN_fastq_raw

    main:
    preprocessing_wf(IN_reference_raw, IN_fastq_raw)

    assembly_wf(IN_fastq_raw)

}