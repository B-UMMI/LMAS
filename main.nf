#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import Helper
import CheckParams

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { preprocessing_wf } from './modules/preprocessing/preprocessing'
include { assembly_wf } from './modules/assembly/assembly'
include { mapping_wf } from './modules/mapping/mapping'
include { postprocessing_wf } from './modules/postprocessing/postprocessing.nf'
include { report_wf } from './modules/report/report.nf'

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

params.version = false
if (params.version){
    println("LMAS version $workflow.manifest.version")
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

Help.start_info(infoMap, "$workflow.start", "$workflow.profile", "$workflow.manifest.version")

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {
    //  Main parameters
    IN_fastq_raw = Channel.fromFilePairs(params.fastq).ifEmpty {exit 1, "No fastq files provided with pattern:'${params.fastq}'"}
    IN_reference_raw = Channel.fromPath(params.reference).ifEmpty {exit 1, "No reference fasta file provided with pattern:'${params.reference}'"}


    // LMAS Steps
    preprocessing_wf(IN_reference_raw, IN_fastq_raw)

    assembly_wf(IN_fastq_raw)

    mapping_wf(assembly_wf.out.all_assemblies, 
               preprocessing_wf.out.triple_reference)

    postprocessing_wf(preprocessing_wf.out.triple_reference,
                      mapping_wf.out.boc_csv, 
                      mapping_wf.out.lx_csv, 
                      mapping_wf.out.nax_csv,
                      mapping_wf.out.ngx_csv,
                      mapping_wf.out.phred_csv,
                      mapping_wf.out.df_csv,
                      mapping_wf.out.paf,
                      mapping_wf.out.misassembly_paf)
 
    report_wf(preprocessing_wf.out.reads_info | collect, 
              mapping_wf.out.stats_global, 
              IN_reference_raw, 
              postprocessing_wf.out.contig_distribution_json, 
              mapping_wf.out.stats_mapping, 
              postprocessing_wf.out.completness_json, 
              postprocessing_wf.out.lx_json, 
              postprocessing_wf.out.phred_json, 
              postprocessing_wf.out.gap_reference_json, 
              postprocessing_wf.out.snp_reference_json, 
              postprocessing_wf.out.gap_boxplot_json, 
              postprocessing_wf.out.misassembly_json, 
              postprocessing_wf.out.misassembly_report_json, 
              postprocessing_wf.out.nax_json, 
              postprocessing_wf.out.ngx_json, 
              postprocessing_wf.out.misassembly_reference_json, 
              postprocessing_wf.out.misassembly_plot_json, 
              assembly_wf.out.all_versions,
              postprocessing_wf.out.snp_report_json)
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
