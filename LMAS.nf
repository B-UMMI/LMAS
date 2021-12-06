#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import Helper
import CheckParams

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { PROCESS_REFERENCE ; PROCESS_READS ; PROCESS_VERSION } from './modules/preprocessing/preprocessing'
include { assembly_wf } from './modules/assembly/assembly'
include { mapping_wf } from './modules/mapping/mapping'
include { COMPILE_REPORT } from './modules/report/report.nf'

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
CollectInitialMetadata.print_metadata(workflow)


process PROCESS_COMPLETNESS {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    path coverage_files 

    output:
    path('*.html') optional true
    path('completness_plots.json'), emit: json

    script:
    template "completness_plot.py"
}

process PLOT_LX {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    path lx_files 
    val(scale) 

    output:
    path('*.html') optional true
    path('lx.json'), emit: json

    script:
    template "lx_plot.py"
}

process PLOT_NAX {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    path nax_files
    val(scale)

    output:
    path('*.html') optional true
    path('nax.json'), emit: json

    script:
    template "nax_plot.py"
}

process PLOT_NGX {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    path ngx_files 
    val(scale) 

    output:
    path('*.html') optional true
    path('ngx.json'), emit: json

    script:
    template "ngx_plot.py"
}

process PROCESS_SHRIMP_PLOT {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    path phred_files 

    output:
    path('*.html') optional true
    path('phred.json'), emit: json

    script:
    template "shrimp_plot.py"
}

process PLOT_CONTIG_DISTRIBUTION {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    path dataframes 

    output:
    path('*.html') optional true
    path('*.json'), emit: json

    script:
    template "plot_contig_size.py"
}

process GAP_ASSESSMENT {

    tag { assembler }
    label 'process_script'
    publishDir "results/$sample_id/stats/"

    input:
    tuple val(sample_id), val(assembler), path(assembly), path(mapping) 
    each reference 

    output:
    path('*_gap_dict.json'), emit: json
    path('*_gaps.csv'), emit: csv

    script:
    template "gap_assessment.py"
}

process PLOT_GAP_BOXPLOT {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    path gap_distance_json 

    output:
    path('*.html') optional true
    path('*gap_distance_histogram.json'), emit: json

    script:
    template "plot_gap_sizes.py"

}

process PLOT_GAP_REFERENCE {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    path gap_coords_dataframes 

    output:
    path('*.html') optional true
    path('*.json'), emit: json

    script:
    template "plot_gap_reference.py"
}

process SNP_ASSESSMENT {

    tag { assembler }
    label 'process_script'

    input:
    tuple val(sample_id), val(assembler), path(assembly), path(mapping) 
    each reference

    output:
    path('*.tsv'), emit: tsv
    path('*_snps.csv'), emit: csv

    script:
    template "snp_assessment.py"
}

process PLOT_SNP_REFERENCE {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    path snp_coords_dataframes 

    output:
    path('*.html') optional true
    path('*.json'), emit: json

    script:
    template "plot_snp.py"
}

process MISASSEMBLY {

    tag { assembler }
    label 'process_script'

    input:
    tuple val(sample_id), val(assembler), path(assembly), path(mapping) 

    output:
    path('*_trace.pkl'), emit: trace_pkl
    path('*_contig_lenght.pkl'), emit: contig_length_pkl
    path('*_misassembly.json'), emit: misassembly_json
    path('*_misassembled_reference.json'), emit: misassembled_reference_json
    path('*_misassembly.csv'), emit: csv

    script:
    template "misassembly.py"

}

process PROCESS_MISASSEMBLY {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    path misassembly_trace 
    path misassembly_contigs
    path report_data 
    path report_per_reference 

    output: 
    path('*.html') optional true
    path('*_misassembly.json'), emit: json
    path('misassembly_report.json'), emit: report_json
    path('misassembly_report_per_ref.json'), emit: reference_json

    script:
    template "process_misassembly.py"

}

process PLOT_MISASSEMBLY {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    path misassembly_dataframes 

    output:
    path('*.html') optional true
    path('*.json'), emit: json

    script:
    template "plot_misassembly.py"

}


/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {

    //  Main parameters
    IN_fastq_raw = Channel.fromFilePairs(params.fastq, size: -1).ifEmpty {exit 1, "No fastq files provided with pattern:'${params.fastq}'"}
    IN_reference_raw = Channel.fromPath(params.reference).ifEmpty {exit 1, "No reference fasta file provided with pattern:'${params.reference}'"}
    IN_MD = Channel.fromPath(params.md).ifEmpty {'No information provided.'}

    //  Optional parameters
    def plot_mode_expected = ['linear', 'log'] as Set
    def plot_parameter_diff = plot_mode_expected - params.plot_scale
    if (plot_parameter_diff.size() > 1) {
        println "[Pipeline warning] Parameter --plot_scale is not valid! Running with default 'log'\n"
        IN_PLOT_SCALE = Channel.from('log')
    } else {
        IN_PLOT_SCALE = Channel.from(params.plot_scale)
    }


    // LMAS Steps
    PROCESS_REFERENCE(IN_reference_raw)
    PROCESS_READS(IN_fastq_raw)

    assembly_wf(IN_fastq_raw)

    PROCESS_VERSION(assembly_wf.out.all_versions)

    mapping_wf(assembly_wf.out.all_assemblies, PROCESS_REFERENCE.out)
 
    PROCESS_COMPLETNESS(mapping_wf.out.boc_csv | collect)
    PLOT_LX(mapping_wf.out.lx_csv, IN_PLOT_SCALE)
    PLOT_NAX(mapping_wf.out.nax_csv, IN_PLOT_SCALE)
    PLOT_NGX(mapping_wf.out.ngx_csv, IN_PLOT_SCALE)
    PROCESS_SHRIMP_PLOT(mapping_wf.out.phred_csv)
    PLOT_CONTIG_DISTRIBUTION(mapping_wf.out.df_csv)

    GAP_ASSESSMENT(mapping_wf.out.paf, PROCESS_REFERENCE.out)
    PLOT_GAP_BOXPLOT(GAP_ASSESSMENT.out.json)
    PLOT_GAP_REFERENCE(GAP_ASSESSMENT.out.csv)
    SNP_ASSESSMENT(mapping_wf.out.paf, PROCESS_REFERENCE.out)
    PLOT_SNP_REFERENCE(SNP_ASSESSMENT.out.csv)
    MISASSEMBLY(mapping_wf.out.paf)
    PROCESS_MISASSEMBLY(MISASSEMBLY.out.trace_pkl | collect, MISASSEMBLY.out.contig_length_pkl | collect, MISASSEMBLY.out.misassembly_json | collect, MISASSEMBLY.out.misassembled_reference_json | collect)
    PLOT_MISASSEMBLY(MISASSEMBLY.out.csv | collect)
    
    pipeline_stats = Channel.fromPath("${workflow.projectDir}/pipeline_stats.txt")
    js = Channel.fromPath("${workflow.projectDir}/resources/main.js.zip")
    lmas_png = Channel.fromPath("${workflow.projectDir}/resources/lmas.zip")
    containers_config = Channel.fromPath("${workflow.projectDir}/configs/containers.config")
    COMPILE_REPORT(PROCESS_READS.out, mapping_wf.out.stats_global, pipeline_stats, js, lmas_png, IN_reference_raw, PLOT_CONTIG_DISTRIBUTION.out.json, mapping_wf.out.stats_mapping, PROCESS_COMPLETNESS.out.json, PLOT_LX.out.json, PROCESS_SHRIMP_PLOT.out.json, PLOT_GAP_REFERENCE.out.json, PLOT_SNP_REFERENCE.out.json, PLOT_GAP_BOXPLOT.out.json, PROCESS_MISASSEMBLY.out.json, PROCESS_MISASSEMBLY.out.report_json, PLOT_NAX.out.json, PLOT_NGX.out.json, PROCESS_VERSION.out, PROCESS_MISASSEMBLY.out.reference_json, PLOT_MISASSEMBLY.out.json, IN_MD, containers_config)
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
