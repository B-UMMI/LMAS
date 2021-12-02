#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import Helper
import CheckParams

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { PROCESS_REFERENCE ; PROCESS_READS ; PROCESS_VERSION } from './modules/processing/processing'
include { assembly_wf } from './modules/assembly/assembly'

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

process FILTER_ASSEMBLY {

    tag { sample_id; assembler }
    label 'process_script'
    publishDir "results/$sample_id/assembly/filtered/"

    input:
    tuple val(sample_id), val(assembler), path(assembly)
    val minLen

    output:
    tuple val(sample_id), val(assembler), path('filtered_*')

    script:
    "reformat.sh in=${assembly} out=filtered_${assembly} minlength=${minLen}"
}

process READ_MAPPING{

    tag { assembler }
    label 'process_mapping'
    publishDir "results/$sample_id/mapping/reads"

    input:
    tuple val(sample_id), val(assembler), path(assembly), path(filtered_assembly)

    output:
    path('*_read_mapping_*.txt') optional true
    tuple val(sample_id), val(assembler), path('*_read_mapping_report.json'), emit: read_mapping_json

    script:
    template "read_mapping.py"
}

process ASSEMBLY_MAPPING{

    tag { sample_id; assembler }
    label 'process_mapping'
    publishDir "results/$sample_id/mapping/assembly"

    input:
    tuple val(sample_id), val(assembler), path(assembly) 
    each path(reference)

    output:
    tuple val(sample_id), val(assembler), path(assembly), file('*.paf')

    script:
    """
    minimap2 --cs -N 50 --secondary=no -t $task.cpus -r 10000 -g 10000 -x asm20 --eqx ${reference} ${assembly} \
    > ${sample_id}_${assembler}.paf
    """

}

//OUT_ASSEMBLY_MAPPING
//.into { IN_ASSEMBLY_MAPPING_FOR_STATS; IN_GAP_ASSESSMENT; IN_SNP_ASSESSMENT; IN_MISASSEMBLY }

// ASSEMBLY STATS GLOBAL
process ASSEMBLY_STATS_GLOBAL {

    tag { assembler }
    label 'process_script'
    publishDir "results/$sample_id/stats/assembly"

    input:
    tuple val(sample_id), val(assembler), path(assembly), path(read_mapping) 

    output:
    path('*report.json'), emit: json
    path('*.csv'), emit: tsv

    script:
    template "assembly_stats_global.py"
}

process PROCESS_ASSEMBLY_STATS_GLOBAL {

    label 'process_script'
    publishDir 'results/stats'

    input:
    path(assembly_stats_global_files)
    path(json_report)

    output:
    path('global_assembly_stats.json')

    script:
    template "process_assembly_stats_global.py"

}

process ASSEMBLY_STATS_MAPPING {

    tag { assembler }
    label 'process_script'
    publishDir "results/$sample_id/stats/"

    input:
    tuple val(sample_id), val(assembler), path(assembly), path(mapping) 
    each reference

    output:
    path('*_report.json'), emit: json
    path('*breadth_of_coverage_contigs.csv'), emit: boc_csv
    path('*_df.csv'), emit: df_csv
    path('*_lx.csv'), emit: lx_csv
    path('*_nax.csv'), emit: nax_csv
    path('*_ngx.csv'), emit: ngx_csv
    path('*_phred.csv'), emit: phred_csv

    script:
    template "assembly_stats_mapping.py"

}

process PROCESS_ASSEMBLY_STATS_MAPPING {

    label 'process_script'
    publishDir 'results/stats/'

    input:
    file json_report

    output:
    path('global_assembly_mapping_stats.json')

    script:
    template "process_assembly_stats_mapping.py"

}

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
    file nax_files
    val(scale)

    output:
    file '*.html' optional true
    file 'nax.json'

    script:
    template "nax_plot.py"
}

process PLOT_NGX {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    file ngx_files 
    val(scale) 

    output:
    file '*.html' optional true
    file 'ngx.json'

    script:
    template "ngx_plot.py"
}

process PROCESS_SHRIMP_PLOT {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    file phred_files 

    output:
    file '*.html' optional true
    file 'phred.json'

    script:
    template "shrimp_plot.py"
}

process PLOT_CONTIG_DISTRIBUTION {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    file dataframes 

    output:
    file '*.html' optional true
    file '*.json'

    script:
    template "plot_contig_size.py"
}

process GAP_ASSESSMENT {

    tag { assembler }
    label 'process_script'
    publishDir "results/$sample_id/stats/"

    input:
    tuple sample_id, assembler, file(assembly), file(mapping) 
    each reference 

    output:
    file '*_gap_dict.json'
    file '*_gaps.csv'

    script:
    template "gap_assessment.py"
}

process PLOT_GAP_BOXPLOT {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    file gap_distance_json 

    output:
    file '*.html' optional true
    file '*gap_distance_histogram.json'

    script:
    template "plot_gap_sizes.py"

}

process PLOT_GAP_REFERENCE {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    file gap_coords_dataframes 

    output:
    file '*.html' optional true
    file '*.json'

    script:
    template "plot_gap_reference.py"
}

process SNP_ASSESSMENT {

    tag { assembler }
    label 'process_script'

    input:
    tuple sample_id, assembler, file(assembly), file(mapping) 
    each reference

    output:
    file '*.tsv'
    file '*_snps.csv'

    script:
    template "snp_assessment.py"
}

process PLOT_SNP_REFERENCE {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    file snp_coords_dataframes 

    output:
    file '*.html' optional true
    file '*.json'

    script:
    template "plot_snp.py"
}

process MISASSEMBLY {

    tag { assembler }
    label 'process_script'

    input:
    tuple sample_id, assembler, file(assembly), file(mapping) 

    output:
    file '*_trace.pkl'
    file '*_contig_lenght.pkl'
    file '*_misassembly.json'
    file '*_misassembled_reference.json'
    file '*_misassembly.csv'

    script:
    template "misassembly.py"

}

process PROCESS_MISASSEMBLY {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    file misassembly_trace 
    file misassembly_contigs
    file report_data 
    file report_per_reference 

    output: 
    file '*.html' optional true
    file '*_misassembly.json'
    file 'misassembly_report.json'
    file 'misassembly_report_per_ref.json'

    script:
    template "process_misassembly.py"

}

process PLOT_MISASSEMBLY {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    file misassembly_dataframes 

    output:
    file '*.html' optional true
    file '*.json'

    script:
    template "plot_misassembly.py"

}

/** Reports
Compiles the reports from every process
**/

process compile_reports {

    label 'process_script'
    publishDir 'report/', mode: 'copy'

    input:
    file reads_json 
    file global_assembly_stats
    file pipeline_stats
    file js
    file lmas_png 
    file reference_file
    file contig_size_distribution 
    file mapping_assembly_stats 
    file completness_plots 
    file lx_plots 
    file shrimp_plots 
    file gap_reference_json 
    file snp_reference_json
    file gap_histogram 
    file plot_misassemblies 
    file misassembly_data 
    file nax_plots 
    file ngx_plots 
    file versions_json 
    file misassembly_per_ref 
    file plot_misassembly_per_ref 
    file about_md 
    file containers_config 

    output:
    file 'pipeline_report*.json'
    file 'index.html'
    file 'main.js'
    file '*.jpg'
    file 'performance_metadata.json'
    file 'reference_metadata.json'

    script:
    template "compile_reports.py"
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
    minLength = Channel.value(params.minLength)

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

    FILTER_ASSEMBLY(assembly_wf.out.all_assemblies, minLength)
    READ_MAPPING(assembly_wf.out.all_assemblies | join(FILTER_ASSEMBLY.out, by:[0,1]))
    ASSEMBLY_MAPPING(FILTER_ASSEMBLY.out, IN_reference_raw)

    ASSEMBLY_STATS_GLOBAL(assembly_wf.out.all_assemblies | join(READ_MAPPING.out.read_mapping_json, by:[0,1]))
    PROCESS_ASSEMBLY_STATS_GLOBAL(ASSEMBLY_STATS_GLOBAL.out.tsv | collect, ASSEMBLY_STATS_GLOBAL.out.json | collect)

    ASSEMBLY_STATS_MAPPING(ASSEMBLY_MAPPING.out, IN_reference_raw)
    PROCESS_ASSEMBLY_STATS_MAPPING(ASSEMBLY_STATS_MAPPING.out.json | collect)

    PROCESS_COMPLETNESS(ASSEMBLY_STATS_MAPPING.out.boc_csv | collect)
    PLOT_LX(ASSEMBLY_STATS_MAPPING.out.lx_csv | collect, IN_PLOT_SCALE)
    PLOT_NAX(ASSEMBLY_STATS_MAPPING.out.nax_csv | collect, IN_PLOT_SCALE)
    PLOT_NGX(ASSEMBLY_STATS_MAPPING.out.ngx_csv | collect, IN_PLOT_SCALE)
    PROCESS_SHRIMP_PLOT(ASSEMBLY_STATS_MAPPING.out.phred_csv | collect)
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
