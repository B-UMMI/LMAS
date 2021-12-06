nextflow.enable.dsl=2

// PROCESSES
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
