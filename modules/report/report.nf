nextflow.enable.dsl=2

// PROCESSES
process COMPILE_VERSIONS {

    label 'process_script'

    input:
    file version

    output:
    file 'versions.json'

    script:
    template "process_versions.py"
}

process COMPILE_REPORT {

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
    file snp_report_json

    output:
    path('pipeline_report*.json')
    path('index.html')
    path('main.js')
    path('*.jpg')
    path('performance_metadata.json')
    path('reference_metadata.json')

    script:
    template "compile_reports.py"
}

workflow report_wf {

    IN_MD = Channel.fromPath(params.md).ifEmpty {'No information provided.'}
    pipeline_stats = Channel.fromPath("${workflow.projectDir}/pipeline_stats.txt")
    js = Channel.fromPath("${workflow.projectDir}/resources/main.js.zip")
    lmas_png = Channel.fromPath("${workflow.projectDir}/resources/lmas.zip")
    containers_config = Channel.fromPath("${workflow.projectDir}/conf/containers.config")
    
    take:
    reads_info
    stats_global
    reference
    plot_contig_distribution
    stats_mapping
    plot_completness
    plot_lx
    plot_phred
    plot_gap_reference
    plot_snp_reference
    plot_gap_boxplot
    misassembly_info
    misassembly_report
    plot_nax
    plot_ngx
    misassembly_reference
    plot_misassembly
    version_file
    snp_report_json

    main:
    COMPILE_VERSIONS(version_file)
    COMPILE_REPORT(reads_info, stats_global, pipeline_stats, js, lmas_png, reference, plot_contig_distribution, stats_mapping, plot_completness, plot_lx, plot_phred, plot_gap_reference, plot_snp_reference, plot_gap_boxplot, misassembly_info, misassembly_report, plot_nax, plot_ngx, COMPILE_VERSIONS.out, misassembly_reference, plot_misassembly, IN_MD, containers_config, snp_report_json)

}