nextflow.enable.dsl=2

// PROCESSES
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

// WORKFLOW
workflow postprocessing_wf {

    def plot_mode_expected = ['linear', 'log'] as Set
    def plot_parameter_diff = plot_mode_expected - params.plot_scale
    if (plot_parameter_diff.size() > 1) {
        println "[Pipeline warning] Parameter --plot_scale is not valid! Running with default 'log'\n"
        IN_PLOT_SCALE = Channel.from('log')
    } else {
        IN_PLOT_SCALE = Channel.from(params.plot_scale)
    }

    take:
    triple_reference
    boc_csv
    lx_csv
    nax_csv
    ngx_csv
    phred_csv
    mapping_df_csv
    paf

    main:
    PROCESS_COMPLETNESS(boc_csv | collect)
    PLOT_LX(lx_csv, IN_PLOT_SCALE)
    PLOT_NAX(nax_csv, IN_PLOT_SCALE)
    PLOT_NGX(ngx_csv, IN_PLOT_SCALE)
    PROCESS_SHRIMP_PLOT(phred_csv)
    PLOT_CONTIG_DISTRIBUTION(mapping_df_csv)
    GAP_ASSESSMENT(paf, triple_reference)
    PLOT_GAP_BOXPLOT(GAP_ASSESSMENT.out.json)
    PLOT_GAP_REFERENCE(GAP_ASSESSMENT.out.csv)
    SNP_ASSESSMENT(paf, triple_reference)
    PLOT_SNP_REFERENCE(SNP_ASSESSMENT.out.csv)
    MISASSEMBLY(paf)
    PROCESS_MISASSEMBLY(MISASSEMBLY.out.trace_pkl | collect, MISASSEMBLY.out.contig_length_pkl | collect, MISASSEMBLY.out.misassembly_json | collect, MISASSEMBLY.out.misassembled_reference_json | collect)
    PLOT_MISASSEMBLY(MISASSEMBLY.out.csv | collect)

    emit:
    completness_json = PROCESS_COMPLETNESS.out.json
    contig_distribution_json = PLOT_CONTIG_DISTRIBUTION.out.json
    lx_json = PLOT_LX.out.json
    phred_json = PROCESS_SHRIMP_PLOT.out.json
    gap_reference_json = PLOT_GAP_REFERENCE.out.json
    snp_reference_json = PLOT_SNP_REFERENCE.out.json
    gap_boxplot_json = PLOT_GAP_BOXPLOT.out.json
    misassembly_json = PROCESS_MISASSEMBLY.out.json
    misassembly_report_json = PROCESS_MISASSEMBLY.out.report_json
    nax_json = PLOT_NAX.out.json
    ngx_json = PLOT_NGX.out.json
    misassembly_reference_json = PROCESS_MISASSEMBLY.out.reference_json
    misassembly_plot_json = PLOT_MISASSEMBLY.out.json
}