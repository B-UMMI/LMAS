process FILTER_ASSEMBLY {

    tag { sample_id; assembler }
    label 'process_script'
    publishDir "results/$sample_id/assembly/filtered/"

    input:
    tuple sample_id, assembler, file(assembly) from TO_FILTER
    val minLen from Channel.value(params.minLength)

    output:
    tuple sample_id, assembler, file('filtered_*') into OUT_FILTERED

    script:
    "reformat.sh in=${assembly} out=filtered_${assembly} minlength=${minLen}"
}

process READ_MAPPING{

    tag { assembler }
    label 'process_mapping'
    publishDir "results/$sample_id/mapping/reads"

    input:
    tuple sample_id, assembler, assembly, filtered_assembly from TO_READ_MAPPING_ALL.join(IN_READ_MAPPING_FILTERED, by: [0,1])

    output:
    file '*_read_mapping_*.txt' optional true
    tuple sample_id, assembler, file('*_read_mapping_report.json') into OUT_READ_MAPPING

    script:
    template "read_mapping.py"
}

process ASSEMBLY_MAPPING{

    tag { sample_id; assembler }
    label 'process_mapping'
    publishDir "results/$sample_id/mapping/assembly"

    input:
    tuple sample_id, assembler, file(assembly) from IN_ASSEMBLY_MAPPING
    each reference from IN_MAPPING_CONTIGS

    output:
    tuple sample_id, assembler, file(assembly), file('*.paf') into OUT_ASSEMBLY_MAPPING

    script:
    """
    minimap2 --cs -N 50 --secondary=no -t $task.cpus -r 10000 -g 10000 -x asm20 --eqx ${reference} ${assembly} \
    > ${sample_id}_${assembler}.paf
    """

}