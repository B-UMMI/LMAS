#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import Helper
import CheckParams

include { LMAS } from './workflows/lmas.nf'
include { LONGLMAS } from './workflows/longlmas.nf'

/*
========================================================================================
    HELP MESSAGE
========================================================================================
*/

// Help message
params.help = false
if (params.help){
    Help.print_help(params)
    exit 0
}

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

Params.check(params)

// Metadata collection for start message
def infoMap = [:]
if (params.containsKey('fastq')) {
    infoMap.put('fastq', file(params.fastq).size())
}
if (params.containsKey('reference')) {
    if (file(params.reference) instanceof LinkedList) {
        infoMap.put('reference', file(params.reference).size())
    } else {
        infoMap.put('fasta', 1)
    }
}

Help.start_info(infoMap, "$workflow.start", "$workflow.profile", "$workflow.manifest.version")

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

//  Main parameters
IN_reference_raw = Channel.fromPath(params.reference).ifEmpty {exit 1, "No reference fasta file provided with pattern:'${params.reference}'"}


 workflow {
    if (params.wf == "default" || params.wf == "Illumina" || params.wf == "illumina") {

        IN_fastq_raw = Channel.fromFilePairs(params.fastq).ifEmpty {exit 1, "No fastq files provided with pattern:'${params.fastq}'"}

        LMAS(IN_reference_raw, IN_fastq_raw)

    } else if (params.wf == "ONT" || params.wf == "ont")  {

        if (params.fastq && !params.ont) {
            IN_fastq_raw = Channel.fromFilePairs(params.fastq).ifEmpty {exit 1, "No fastq files provided with pattern:'${params.fastq}'"}
        } else if (!params.fastq && params.ont) {
            IN_fastq_raw = Channel.fromFilePairs(params.ont).ifEmpty {exit 1, "No fastq files provided with pattern:'${params.ont}'"}
        } else {
            exit 1, "Input reads should be provided with --fastq or --ont"
        }

        LONGLMAS(IN_reference_raw, IN_fastq_raw)

    } else if (params.wf == "Hybrid" || params.wf == "hybrid")  {

        IN_fastq_raw = Channel.fromFilePairs(params.fastq).ifEmpty {exit 1, "No fastq files provided with pattern:'${params.fastq}'"}
        IN_ont_raw = Channel.fromFilePairs(params.ont).ifEmpty {exit 1, "No fastq files provided with pattern:'${params.ont}'"}

        HYBRIDLMAS(IN_reference_raw, IN_fastq_raw, IN_ont_raw)

    } else {
        exit 1, "Unrecogized --wf parameter: '${params.wf}'"
    }
}

workflow.onComplete {
  // Display complete message
  log.info "Completed at: " + workflow.complete
  log.info "Duration    : " + workflow.duration
  log.info "Success     : " + workflow.success
  log.info "Exit status : " + workflow.exitStatus
}

workflow.onError {
  // Display error message
  log.info "Workflow execution stopped with the following message:"
  log.info "  " + workflow.errorMessage
}
