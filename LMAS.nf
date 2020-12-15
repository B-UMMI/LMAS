#!/usr/bin/env nextflow

import Helper
import CollectInitialMetadata

// Pipeline version
if (workflow.commitId){
    version = "0.1 $workflow.revision"
} else {
    version = "0.1 (local version)"
}

params.help = false
if (params.help){
    Help.print_help(params)
    exit 0
}

def infoMap = [:]
if (params.containsKey("fastq")){
    infoMap.put("fastq", file(params.fastq).size())
}
if (params.containsKey("reference")){
    if (file(params.reference) instanceof LinkedList){
        infoMap.put("reference", file(params.reference).size())
    } else {
        infoMap.put("fasta", 1)
    }
}

Help.start_info(infoMap, "$workflow.start", "$workflow.profile", version)
CollectInitialMetadata.print_metadata(workflow)

/*
Workflow Start!
*/

// MAIN PARAMETERS
//      FastQ
if (params.fastq instanceof Boolean){
    exit 1, "'fastq' must be a path pattern. Provided value:'$params.fastq'"
    }
if (!params.fastq){ exit 1, "'fastq' parameter missing"}
// size: -1 -> allows for single and paired-end files to be passed through. Change if necessary
IN_fastq_raw = Channel.fromFilePairs(params.fastq, size: -1).ifEmpty {
    exit 1, "No fastq files provided with pattern:'${params.fastq}'" }

//      Reference
if (params.reference instanceof Boolean){
    exit 1, "'reference' must be a path pattern. Provided value: '$params.reference'"
}
if (!params.reference){ exit 1, "'reference' parameter missing"}
IN_reference_raw = Channel.fromPath(params.reference).ifEmpty {
    exit 1, "No reference fasta file provided with pattern:'${params.reference}'" }


// SET CHANNELS FOR ASSEMBLERS
IN_fastq_raw.into{
    IN_BCALM2;
    IN_GATB_MINIA_PIPELINE;
    IN_MINIA;
    IN_MEGAHIT;
    IN_METASPADES;
    IN_UNICYCLER;
    IN_IDBA;
    IN_SPADES;
    IN_SKESA;
    IN_VELVETOPTIMIZER;
    IN_PANDASEQ;
    IN_TO_MAP} //mapping channel - minimap2

// SET CHANNELS FOR REFERENCE
IN_reference_raw.into{IN_MAPPING_CONTIGS; IN_ASSEMBLY_STATS_MAPPING}

// ASSEMBLERS
//      BCALM 2
if ( !params.bcalmKmerSize.toString().isNumber() ){
    exit 1, "'bcalmKmerSize' parameter must be a number. Provided value: '${params.bcalmKmerSize}'"
}

process BCALM2 {
    tag {sample_id}
    publishDir "results/assembly/bcalm2/"

    input:
    set sample_id, file(fastq) from IN_BCALM2
    val KmerSize from Channel.value(params.bcalmKmerSize)

    output:
    set sample_id, val("BCALM2"), file("*_BCALM2.fasta") into OUT_BCALM2

    script:
    """
    ls -1 $fastq  > list_reads
    {
        bcalm -in list_reads -out ${sample_id} -kmer-size $KmerSize
        mv ${sample_id}.unitigs.fa  ${sample_id}_BCALM2.fasta
        echo pass > .status
    } || {
        echo fail > .status
    }
    # remove temp files
    rm list_reads *.fa || true
    """
}


//      GATB MINIA Pipeline
IN_GATB_kmers = Channel.value(params.gatbkmer)
IN_GATB_besst_iter = Channel.value(params.gatb_besst_iter)
GATB_error_correction = params.GATB_error_correction ? "true" : "false"
IN_error_correction = Channel.value(GATB_error_correction)

process GATBMINIAPIPELINE {
    tag {sample_id}
    publishDir 'results/assembly/GATBMiniaPipeline/'

    input:
    set sample_id, file(fastq_pair) from IN_GATB_MINIA_PIPELINE
    val kmer_list from IN_GATB_kmers
    val do_error_correction from GATB_error_correction
    val besst_iter from IN_GATB_besst_iter

    output:
    set sample_id, val("GATBMiniaPipeline"), file('*_GATBMiniaPipeline.fasta') into OUT_GATB

    script:
    """
    {
        if [ $do_error_correction ];
        then
            gatb -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} --kmer-sizes ${kmer_list} -o ${sample_id}_GATBMiniaPipeline
        else
            gatb -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} --kmer-sizes ${kmer_list} -o ${sample_id}_GATBMiniaPipeline --no-error-correction
        fi

        link=\$(readlink ${sample_id}_GATBMiniaPipeline.fasta) && rm ${sample_id}_GATBMiniaPipeline.fasta && mv \$link ${sample_id}_GATBMiniaPipeline.fasta

        # rm temp dirs

        echo pass > .status
    } || {
        echo fail > .status
    }
    rm -r *_GATBMiniaPipeline.lib* *_GATBMiniaPipeline_besst *.unitigs* *contigs.fa *.h5 || true
    rm *list_reads* || true
    """
}


//      MINIA
IN_MINIA_kmer = Channel.value(params.miniakmer)

process MINIA {
    tag {sample_id}
    publishDir 'results/assembly/MINIA/'

    input:
    set sample_id, file(fastq) from IN_MINIA
    val kmer from IN_MINIA_kmer

    output:
    set sample_id, val("MINIA"), file('*_minia.fasta') into OUT_MINIA

    script:
    """
    {
        ls -1 $fastq  > list_reads
        minia -in list_reads -out ${sample_id}_minia.fasta -nb-cores $task.cpu
        mv ${sample_id}_minia.fasta.contigs.fa ${sample_id}_minia.fasta
        echo pass > .status
    } || {
        echo fail > .status
    }
    rm list_reads *.unitigs.* *.h5 || true
    """
}


//      MEGAHIT
IN_megahit_kmers = Channel.value(params.megahitKmers)

process MEGAHIT {
    tag { sample_id }
    publishDir 'results/assembly/MEGAHIT/', pattern: '*_megahit*.fasta'

    input:
    set sample_id, file(fastq_pair) from IN_MEGAHIT
    val kmers from IN_megahit_kmers

    output:
    set sample_id, val("MEGAHIT"), file('*_MEGAHIT.fasta') into OUT_MEGAHIT

    script:
    """
    {
        /NGStools/megahit/bin/megahit --num-cpu-threads $task.cpus -o megahit --k-list $kmers -1 ${fastq_pair[0]} -2 ${fastq_pair[1]}
        mv megahit/final.contigs.fa ${sample_id}_MEGAHIT.fasta
        echo pass > .status
    } || {
        echo fail > .status
    }
    rm -r megahit || true
    """

}


//      METASPADES
if ( params.metaspadesKmers.toString().split(" ").size() <= 1 ){
    if (params.metaspadesKmers.toString() != 'auto'){
        exit 1, "'metaspadesKmers' parameter must be a sequence of space separated numbers or 'auto'. Provided value: ${params.metaspadesKmers}"
    }
}
IN_metaspades_kmers = Channel.value(params.metaspadesKmers)

process METASPADES {
    tag { sample_id }
    publishDir 'results/assembly/metaSPAdes/'

    input:
    set sample_id, file(fastq_pair) from IN_METASPADES
    val kmers from IN_metaspades_kmers

    output:
    set sample_id, val("metaSPAdes"), file('*_metaspades.fasta') into OUT_METASPADES

    script:
    """
    {
        metaspades.py --only-assembler --threads $task.cpus -k $kmers -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} -o metaspades
        mv metaspades/contigs.fasta ${sample_id}_metaspades.fasta
        echo pass > .status
    } || {
        echo fail > .status
    }
    rm -r metaspades || true
    """
}


//      UNICYCLER
process UNICYCLER {
    tag { sample_id }
    publishDir 'results/assembly/unicycler'

    input:
    set sample_id, file(fastq_pair) from IN_UNICYCLER

    output:
    set sample_id, val("Unicycler"), file('*_unicycler.fasta') into OUT_UNICYCLER

    script:
    """
    {
        unicycler -t $task.cpus -o . --no_correct --no_pilon -1 ${fastq_pair[0]} -2 ${fastq_pair[1]}
        mv assembly.fasta ${sample_id}_unicycler.fasta
        echo pass > .status
    } || {
        echo fail > .status
    }
    rm *best_spades_graph* *overlaps_removed* *bridges_applied* *final_clean* || true
    """
}


//      SPADES
if ( params.spadesKmers.toString().split(" ").size() <= 1 ){
    if (params.spadesKmers.toString() != 'auto'){
        exit 1, "'spadesKmers' parameter must be a sequence of space separated numbers or 'auto'. Provided value: ${params.spadesKmers}"
    }
}
IN_spades_kmers = Channel.value(params.spadesKmers)

process SPADES {
    tag { sample_id }
    publishDir 'results/assembly/SPAdes/', pattern: '*.fasta'

    input:
    set sample_id, file(fastq_pair) from IN_SPADES
    val kmers from IN_spades_kmers

    output:
    set sample_id, val("SPAdes"), file('*_spades.fasta') into OUT_SPADES

    script:
    """
    {
        spades.py --only-assembler --threads $task.cpus -k $kmers -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} -o spades
        mv spades/contigs.fasta ${sample_id}_spades.fasta
    } || {
        echo fail > .status
    }
    rm -r spades || true
    """
}


//      SKESA
process SKESA {
    tag { sample_id }
    publishDir 'results/assembly/SKESA/'

    input:
    set sample_id, file(fastq_pair) from IN_SKESA

    output:
    set sample_id, val("SKESA"), file('*_skesa.fasta') into OUT_SKESA

    script:
    """
    {
        skesa --cores $task.cpus --memory $task.memory --use_paired_ends --contigs_out ${sample_id}_skesa.fasta --fastq ${fastq_pair[0]} ${fastq_pair[1]}
        echo pass > .status
    } || {
        echo fail > .status
    }
    """
}


//      PANDASEQ
process PANDASEQ {
    tag { sample_id }
    publishDir 'results/assembly/pandaseq'

    input:
    set sample_id, file(fastq_pair) from IN_PANDASEQ

    output:
    set sample_id, val("Pandaseq"), file('*pandaseq.fasta') into OUT_PANDASEQ

    script:
    """
    cp -r /NGStools/pandaseq pandaseq/
    {
        ./pandaseq/pandaseq -T $task.cpus -w ${sample_id}_pandaseq.fasta -f ${fastq_pair[0]} -r ${fastq_pair[1]} -B
        echo pass > .status
    } || {
        echo fail > .status
    }
    rm -r pandaseq || true
    """
}


//      VELVETOPTIMIZER
process VELVETOPTIMIZER {
    tag { sample_id }
    publishDir 'results/assembly/VelvetOtimiser'

    input:
    set sample_id, file(fastq_pair) from IN_VELVETOPTIMIZER

    output:
    set sample_id, val("VelvetOptimizer"), file('*.fasta') into OUT_VELVETOPTIMIZER

    script:
    """
    {
        VelvetOptimiser.pl -v -s $params.velvetoptimizer_hashs -e $params.velvetoptimizer_hashe -t $task.cpus \
        -f '-shortPaired -fastq.gz -separate ${fastq_pair[0]} ${fastq_pair[1]}'
        mv auto_data*/contigs.fa ${sample_id}_velvetoptimizer.fasta
        echo pass > .status
    } || {
        echo fail > .status
    }
    rm -r auto_data* || true
    """
}


//      IDBA
process reformat_IDBA {
    tag { sample_id }

    input:
    set sample_id, file(fastq_pair) from IN_IDBA

    output:
    set sample_id, file('*.fasta') into REFORMAT_IDBA

    script:
    "reformat.sh in=${fastq_pair[0]} in2=${fastq_pair[1]} out=${sample_id}_reads.fasta"
}

process IDBA {
    tag { sample_id }
    publishDir 'results/assembly/IDBA-UD/'

    input:
    set sample_id, file(fasta_reads_single) from  REFORMAT_IDBA

    output:
    set sample_id, val("IDBA-UD"), file('*_IDBA-UD.fasta') into OUT_IDBA

    script:
    """
    {
        idba_ud -l ${fasta_reads_single} --num_threads $task.cpus -o .
        mv contig.fa ${sample_id}_IDBA-UD.fasta
        echo pass > .status
    } || {
        echo fail > .status
    }
    rm begin align-* contig-* graph-* kmer local-* || true
    """
}

// ASSEMBLY COLLECTION
OUT_BCALM2.mix(OUT_GATB,
                  OUT_MINIA,
                  OUT_MEGAHIT,
                  OUT_METASPADES,
                  OUT_UNICYCLER,
                  OUT_SPADES,
                  OUT_SKESA,
                  OUT_PANDASEQ,
                  OUT_VELVETOPTIMIZER,
                  OUT_IDBA).set{ALL_ASSEMBLERS}

ALL_ASSEMBLERS.into{ TO_FILTER; TO_GLOBAL_STATS; TO_READ_MAPPING}

THRESHOLD = Channel.value(params.mapped_reads_threshold)
// READ MAPPING
process READ_MAPPING{

    tag { assembler }

    publishDir 'results/stats/${sample_id}/'

    input:
    set sample_id, assembler, assembly from TO_READ_MAPPING
    each THRESHOLD

    output:
    file("*_read_mapping.txt")
    set sample_id, assembler, file("*_read_mapping_report.json") into OUT_READ_MAPPING

    script:
    template "read_mapping.py"
}

// ASSEMBLY STATS GLOBAL
process ASSEMBLY_STATS_GLOBAL {
    tag { assembler }

    publishDir 'results/stats/assembly/'

    input:
    set sample_id, assembler, file(assembly), file(read_mapping) from TO_GLOBAL_STATS.join(OUT_READ_MAPPING, by: [0,1])

    output:
    file "*report.json" into OUT_ASSEMBLY_STATS_GLOBAL_JSON
    file "*.csv" into OUT_ASSEMBLY_STATS_GLOBAL_TSV

    script:
    template "assembly_stats_global.py"
}


process PROCESS_ASSEMBLY_STATS_GLOBAL {

    publishDir 'results/stats/'

    input:
    file assembly_stats_global_files from OUT_ASSEMBLY_STATS_GLOBAL_TSV.collect()
    file json_report from OUT_ASSEMBLY_STATS_GLOBAL_JSON.collect()

    output:
    file "*.csv"
    file "global_assembly_stats.json" into PROCESS_ASSEMBLY_STATS_GLOBAL_OUT

    script:
    template "process_assembly_stats_global.py"

}

// FILTER ASSEMBLY
IN_minLen = Channel.value(params.minLength)

process FILTER_ASSEMBLY {

    tag {sample_id; assembler}
    publishDir 'results/assembly/filtered/'

    input:
    set sample_id, assembler, file(assembly) from TO_FILTER
    val minLen from IN_minLen

    output:
    set sample_id, assembler, file('filtered_*') into OUT_FILTERED

    script:
    "reformat.sh in=${assembly} out=filtered_${assembly} minlength=${minLen}"
}

// ASSEMBLY MAPPING
process ASSEMBLY_MAPPING{

    tag { sample_id; assembler }

    publishDir 'results/mapping/${sample_id}/'

    input:
    set sample_id, assembler, file(assembly) from OUT_FILTERED
    each reference from IN_MAPPING_CONTIGS

    output:
    set sample_id, assembler, file(assembly), file("*.paf") into OUT_ASSEMBLY_MAPPING

    script:
    "minimap2 --cs -t $task.cpus -r 10000 -g 10000 -x asm20 --eqx ${reference} ${assembly} > ${sample_id}_${assembler}.paf"

}

OUT_ASSEMBLY_MAPPING.into{ IN_ASSEMBLY_STATS_MAPPING; IN_GAP_ASSESSMENT}

process ASSEMBLY_STATS_MAPPING {

    tag { assembler }

    publishDir 'results/stats/'

    input:
    set sample_id, assembler, file(assembly), file(mapping) from IN_ASSEMBLY_STATS_MAPPING
    each reference from IN_ASSEMBLY_STATS_MAPPING

    output:
    file("*_report.json") into OUT_ASSEMBLY_STATS_MAPPING_JSON
    file("*breadth_of_coverage_contigs.csv") into OUT_COVERAGE_PER_CONTIG
    file "*_df.csv" into OUT_DF_ASSEMBLY_STATS_MAPPING
    file("*_c90.csv") into OUT_C90
    file("*_phred.csv") into OUT_PHRED

    script:
    template "assembly_stats_mapping.py"

}

process PROCESS_ASSEMBLY_STATS_MAPPING {

    publishDir 'results/stats/'

    input:
    file json_report from OUT_ASSEMBLY_STATS_MAPPING_JSON.collect()

    output:
    file "global_assembly_mapping_stats.json" into PROCESS_ASSEMBLY_STATS_MAPPING_OUT

    script:
    template "process_assembly_stats_mapping.py"

}

process PROCESS_COMPLETNESS {

    publishDir 'results/plots/', pattern: "*.html"

    input:
    file coverage_files from OUT_COVERAGE_PER_CONTIG.collect()

    output:
    file("*.html")
    file("completness_plots.json") into PLOT_PROCESS_COMPLETNESS

    script:
    template "completness_plot.py"
}

process PROCESS_C90 {

    publishDir 'results/plots/', pattern: "*.html"

    input:
    file c90_files from OUT_C90.collect()

    output:
    file("*.html")
    file("c90.json") into PLOT_C90

    script:
    template "c90_plot.py"
}

process PROCESS_SHRIMP_PLOT {

    publishDir 'results/plots/', pattern: "*.html"

    input:
    file phred_files from OUT_PHRED.collect()

    output:
    file("*.html")
    file("phred.json") into PLOT_PHRED

    script:
    template "shrimp_plot.py"
}


process PLOT_CONTIG_DISTRIBUTION {

    publishDir 'results/plots/', pattern: "*.html"

    input:
    file dataframes from OUT_DF_ASSEMBLY_STATS_MAPPING.collect()

    output:
    file("*.html")
    file("*.json") into PLOT_CONTIG_DISTRIBUTION

    script:
    template "plot_contig_size.py"
}

process GAP_ASSESSMENT {

    input:
    set sample_id, assembler, file(assembly), file(mapping) from IN_GAP_ASSESSMENT

    output:
    file("*_gap_dict.json") into OUT_GAP_DISTANCE

    script:
    template "gap_assessment.py"
}

process PLOT_GAP_HISTOGRAM {

    publishDir 'results/plots/', pattern: "*.html"

    input:
    file gap_distance_json from OUT_GAP_DISTANCE.collect()

    output:
    file("*.html")
    file("gap_distance_histogram.json") into OUT_GAP_HISTOGRAM

    script:
    template "plot_gap_distance.py"

}

/** Reports
Compiles the reports from every process
**/

// WORK IN PROGRESS!

OUT_ASSEMBLY_STATS_GLOBAL_JSON.set{master_report}

process compile_reports {

    publishDir "report/", mode: "copy"

    input:
    file global_assembly_stats from PROCESS_ASSEMBLY_STATS_GLOBAL_OUT
    file pipeline_stats from Channel.fromPath("${workflow.projectDir}/pipeline_stats.txt")
    file js from Channel.fromPath("${workflow.projectDir}/resources/main.js.zip")
    file lmas_png from Channel.fromPath("${workflow.projectDir}/resources/lmas.zip")
    file reference_file from Channel.fromPath(params.reference)
    file contig_size_distribution from PLOT_CONTIG_DISTRIBUTION
    file mapping_assembly_stats from PROCESS_ASSEMBLY_STATS_MAPPING_OUT
    file completness_plots from PLOT_PROCESS_COMPLETNESS
    file c90_plots from PLOT_C90
    file shrimp_plots from PLOT_PHRED

    output:
    file "pipeline_report.json"
    file "index.html"
    file "main.js"
    file "lmas.png"

    script:
    template "compile_reports.py"
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