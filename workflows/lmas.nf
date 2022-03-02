/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { preprocessing_wf } from '../modules/preprocessing/preprocessing'
include { assembly_wf } from '../modules/assembly/assembly'
include { mapping_wf } from '../modules/mapping/mapping'
include { postprocessing_wf } from '../modules/postprocessing/postprocessing.nf'
include { report_wf } from '../modules/report/report.nf'

/*
========================================================================================
    RUN WORKFLOW
========================================================================================
*/

workflow LMAS {

    take:
    IN_reference_raw
    IN_fastq_raw

    main:
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
                      mapping_wf.out.paf)

    report_wf(preprocessing_wf.out.reads_info, 
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
            assembly_wf.out.all_versions)
}