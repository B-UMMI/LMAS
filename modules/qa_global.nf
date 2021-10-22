process ASSEMBLY_STATS_GLOBAL {

    tag { assembler }
    label 'process_script'
    publishDir "results/$sample_id/stats/assembly"

    input:
    tuple sample_id, assembler, file(assembly), file(read_mapping) from TO_GLOBAL_STATS.join(OUT_READ_MAPPING, by: [0,1])

    output:
    file '*report.json' into OUT_ASSEMBLY_STATS_GLOBAL_JSON
    file '*.csv' into OUT_ASSEMBLY_STATS_GLOBAL_TSV

    script:
    template "assembly_stats_global.py"
}

process PROCESS_ASSEMBLY_STATS_GLOBAL {

    label 'process_script'
    publishDir 'results/stats'

    input:
    file assembly_stats_global_files from OUT_ASSEMBLY_STATS_GLOBAL_TSV.collect()
    file json_report from OUT_ASSEMBLY_STATS_GLOBAL_JSON.collect()

    output:
    file 'global_assembly_stats.json' into PROCESS_ASSEMBLY_STATS_GLOBAL_OUT

    script:
    template "process_assembly_stats_global.py"

}

process PLOT_CONTIG_DISTRIBUTION {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    file dataframes from OUT_DF_ASSEMBLY_STATS_MAPPING.collect()

    output:
    file '*.html' optional true
    file '*.json' into PLOT_CONTIG_DISTRIBUTION

    script:
    template "plot_contig_size.py"
}

process PLOT_GAP_BOXPLOT {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    file gap_distance_json from OUT_GAP_DISTANCE.collect()

    output:
    file '*.html' optional true
    file '*gap_distance_histogram.json' into OUT_GAP_HISTOGRAM

    script:
    template "plot_gap_sizes.py"

}

