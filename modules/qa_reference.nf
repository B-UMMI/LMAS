process ASSEMBLY_STATS_MAPPING {

    tag { assembler }
    label 'process_script'
    publishDir "results/$sample_id/stats/"

    input:
    tuple sample_id, assembler, file(assembly), file(mapping) from IN_ASSEMBLY_MAPPING_FOR_STATS
    each reference from IN_ASSEMBLY_STATS_MAPPING

    output:
    file '*_report.json' into OUT_ASSEMBLY_STATS_MAPPING_JSON
    file '*breadth_of_coverage_contigs.csv' into OUT_COVERAGE_PER_CONTIG
    file '*_df.csv' into OUT_DF_ASSEMBLY_STATS_MAPPING
    file '*_lx.csv' into OUT_LX_PLOT
    file '*_nax.csv' into OUT_NAX_PLOT
    file '*_ngx.csv' into OUT_NGX_PLOT
    file '*_phred.csv' into OUT_PHRED

    script:
    template "assembly_stats_mapping.py"

}

process PROCESS_ASSEMBLY_STATS_MAPPING {

    label 'process_script'
    publishDir 'results/stats/'

    input:
    file json_report from OUT_ASSEMBLY_STATS_MAPPING_JSON.collect()

    output:
    file 'global_assembly_mapping_stats.json' into PROCESS_ASSEMBLY_STATS_MAPPING_OUT

    script:
    template "process_assembly_stats_mapping.py"

}

process PROCESS_COMPLETNESS {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    file coverage_files from OUT_COVERAGE_PER_CONTIG.collect()

    output:
    file '*.html' optional true
    file 'completness_plots.json' into PLOT_PROCESS_COMPLETNESS

    script:
    template "completness_plot.py"
}

process PLOT_LX {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    file lx_files from OUT_LX_PLOT.collect()
    val(scale) from IN_PLOT_SCALE_1

    output:
    file '*.html' optional true
    file 'lx.json' into PLOT_LX

    script:
    template "lx_plot.py"
}

process PLOT_NAX {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    file nax_files from OUT_NAX_PLOT.collect()
    val(scale) from IN_PLOT_SCALE_2

    output:
    file '*.html' optional true
    file 'nax.json' into PLOT_NAX

    script:
    template "nax_plot.py"
}

process PLOT_NGX {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    file ngx_files from OUT_NGX_PLOT.collect()
    val(scale) from IN_PLOT_SCALE_3

    output:
    file '*.html' optional true
    file 'ngx.json' into PLOT_NGX

    script:
    template "ngx_plot.py"
}

process PROCESS_SHRIMP_PLOT {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    file phred_files from OUT_PHRED.collect()

    output:
    file '*.html' optional true
    file 'phred.json' into PLOT_PHRED

    script:
    template "shrimp_plot.py"
}


process GAP_ASSESSMENT {

    tag { assembler }
    label 'process_script'
    publishDir "results/$sample_id/stats/"

    input:
    tuple sample_id, assembler, file(assembly), file(mapping) from IN_GAP_ASSESSMENT
    each reference from IN_GAP_STATS

    output:
    file '*_gap_dict.json' into OUT_GAP_DISTANCE
    file '*_gaps.csv' into OUT_GAP_PLOT_REF

    script:
    template "gap_assessment.py"
}


process PLOT_GAP_REFERENCE {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    file gap_coords_dataframes from OUT_GAP_PLOT_REF.collect()

    output:
    file '*.html' optional true
    file '*.json' into OUT_GAP_REFERENCE

    script:
    template "plot_gap_reference.py"
}

process SNP_ASSESSMENT {

    tag { assembler }
    label 'process_script'

    input:
    tuple sample_id, assembler, file(assembly), file(mapping) from IN_SNP_ASSESSMENT
    each reference from IN_SNP_STATS

    output:
    file '*.tsv'
    file '*_snps.csv' into OUT_SNP_PLOT_REF

    script:
    template "snp_assessment.py"
}

process PLOT_SNP_REFERENCE {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    file snp_coords_dataframes from OUT_SNP_PLOT_REF.collect()

    output:
    file '*.html' optional true
    file '*.json' into OUT_SNP_REFERENCE

    script:
    template "plot_snp.py"
}

process MISASSEMBLY {

    tag { assembler }
    label 'process_script'

    input:
    tuple sample_id, assembler, file(assembly), file(mapping) from IN_MISASSEMBLY

    output:
    file '*_trace.pkl' into OUT_MISASSEMBLY_TRACE
    file '*_contig_lenght.pkl' into OUT_MISASSEMBLY_CONTIGS
    file '*_misassembly.json' into MISASSEMBLY_REPORT
    file '*_misassembled_reference.json' into MISASSEMBLY_DICTIONARY
    file '*_misassembly.csv' into PLOT_MISASSEMBLY_REF

    script:
    template "misassembly.py"

}

process PROCESS_MISASSEMBLY {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    file misassembly_trace from OUT_MISASSEMBLY_TRACE.collect()
    file misassembly_contigs from OUT_MISASSEMBLY_CONTIGS.collect()
    file report_data from MISASSEMBLY_REPORT.collect()
    file report_per_reference from MISASSEMBLY_DICTIONARY.collect()

    output: 
    file '*.html' optional true
    file '*_misassembly.json' into OUT_MISASSEMBLY_PLOT
    file 'misassembly_report.json' into OUT_MISASSEMBLY_REPORT
    file 'misassembly_report_per_ref.json' into MISASSEMBLY_PER_REF

    script:
    template "process_misassembly.py"

}

process PLOT_MISASSEMBLY {

    label 'process_script'
    publishDir 'results/plots/', pattern: '*.html'

    input:
    file misassembly_dataframes from PLOT_MISASSEMBLY_REF.collect()

    output:
    file '*.html' optional true
    file '*.json' into OUT_MISASSEMBLY_REFERENCE

    script:
    template "plot_misassembly.py"

}