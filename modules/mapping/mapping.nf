nextflow.enable.dsl=2

// PROCESSES
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
    tuple val(sample_id), val(assembler), path(assembly), path('*.paf')

    script:
    """
    minimap2 --cs -N 50 --secondary=yes -t $task.cpus -r 10000 -g 10000 -x asm20 --eqx ${reference} ${assembly} \
    > ${sample_id}_${assembler}.paf
    """
}

process ASSEMBLY_MAPPING_MISASSEMBLY{

    tag { sample_id; assembler }
    label 'process_mapping'
    publishDir "results/$sample_id/mapping/assembly"

    input:
    tuple val(sample_id), val(assembler), path(assembly) 
    each path(reference)

    output:
    tuple val(sample_id), val(assembler), path(assembly), path('*.paf')

    script:
    """
    minimap2 --cs -N 50 --secondary=no -t $task.cpus -r 10000 -g 10000 -x asm20 --eqx ${reference} ${assembly} \
    > ${sample_id}_${assembler}.paf
    """
}

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


// WORKFLOW
workflow mapping_wf {   

    minLength = Channel.value(params.minLength) 

    take:
    assembly
    triple_reference

    main:
    FILTER_ASSEMBLY(assembly, minLength)
    READ_MAPPING(assembly | join(FILTER_ASSEMBLY.out, by:[0,1]))
    ASSEMBLY_MAPPING(FILTER_ASSEMBLY.out, triple_reference)
    ASSEMBLY_MAPPING_MISASSEMBLY(FILTER_ASSEMBLY.out, triple_reference)
    ASSEMBLY_STATS_GLOBAL(assembly | join(READ_MAPPING.out.read_mapping_json, by:[0,1]))
    PROCESS_ASSEMBLY_STATS_GLOBAL(ASSEMBLY_STATS_GLOBAL.out.tsv | collect, ASSEMBLY_STATS_GLOBAL.out.json | collect)
    ASSEMBLY_STATS_MAPPING(ASSEMBLY_MAPPING.out, triple_reference)
    PROCESS_ASSEMBLY_STATS_MAPPING(ASSEMBLY_STATS_MAPPING.out.json | collect)

    emit:
    paf = ASSEMBLY_MAPPING.out
    misassembly_paf = ASSEMBLY_MAPPING_MISASSEMBLY.out
    stats_global = PROCESS_ASSEMBLY_STATS_GLOBAL.out
    stats_mapping = PROCESS_ASSEMBLY_STATS_MAPPING.out
    boc_csv = ASSEMBLY_STATS_MAPPING.out.boc_csv | collect
    lx_csv = ASSEMBLY_STATS_MAPPING.out.lx_csv | collect
    nax_csv = ASSEMBLY_STATS_MAPPING.out.nax_csv | collect
    ngx_csv = ASSEMBLY_STATS_MAPPING.out.ngx_csv | collect
    phred_csv = ASSEMBLY_STATS_MAPPING.out.phred_csv | collect
    df_csv = ASSEMBLY_STATS_MAPPING.out.df_csv | collect
}