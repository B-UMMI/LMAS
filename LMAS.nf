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
    IN_PROCESS_READS;
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
    IN_TO_MAP} //mapping channel - minimap2

// TRIPLE THE REFERENCE REPLICONS
process PROCESS_REFERENCE{
    input:
    file reference_fasta from IN_reference_raw

    output:
    file("triple_reference.fasta") into OUT_REFERENCE_TRIPLE

    script:
    template "process_reference.py"
}

// SET CHANNELS FOR REFERENCE
OUT_REFERENCE_TRIPLE.into{IN_MAPPING_CONTIGS; IN_ASSEMBLY_STATS_MAPPING; IN_GAP_STATS; IN_SNP_STATS; COMPILE_REPORTS_REF}


// ASSEMBLERS
//      BCALM 2
if ( !params.bcalmKmerSize.toString().isNumber() ){
    exit 1, "'bcalmKmerSize' parameter must be a number. Provided value: '${params.bcalmKmerSize}'"
}

process PROCESS_READS{
    tag {sample_id}

    input:
    tuple sample_id, file(fastq) from IN_PROCESS_READS

    output:
    file("*_reads_report.json") into PROCESS_READS

    script:
    template "process_reads.py"
}

process BCALM2 {
    tag {sample_id}
    publishDir "results/$sample_id/assembly/bcalm2/"

    input:
    tuple sample_id, file(fastq) from IN_BCALM2
    val KmerSize from Channel.value(params.bcalmKmerSize)

    output:
    tuple sample_id, val("BCALM2"), file("*_BCALM2.fasta") into OUT_BCALM2
    file(".*version") into BCALM2_VERSION

    script:
    """
    ls -1 $fastq  > list_reads
    bcalm -version | head -n 1 | awk -F ', ' '{print \$2}' | awk -F ' ' '{print \$2}' | awk -F 'v' '{print \$2}' > .${sample_id}_BCALM2_version
    {
        bcalm -in list_reads -out ${sample_id} -kmer-size $KmerSize
        mv ${sample_id}.unitigs.fa  ${sample_id}_BCALM2.fasta
        echo pass > .status
    } || {
        echo fail > .status
        :> ${sample_id}_BCALM2.fasta
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
    publishDir "results/$sample_id/assembly/GATBMiniaPipeline/"

    input:
    tuple sample_id, file(fastq_pair) from IN_GATB_MINIA_PIPELINE
    val kmer_list from IN_GATB_kmers
    val do_error_correction from GATB_error_correction
    val besst_iter from IN_GATB_besst_iter

    output:
    tuple sample_id, val("GATBMiniaPipeline"), file('*_GATBMiniaPipeline.fasta') into OUT_GATB
    file(".*version") into GATB_VERSION

    script:
    """
    echo '' > .${sample_id}_GATBMiniaPipeline_version
    {
        if [ $do_error_correction ];
        then
            gatb -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} --kmer-sizes ${kmer_list} -o ${sample_id}_GATBMiniaPipeline --no-scaffolding
        else
            gatb -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} --kmer-sizes ${kmer_list} -o ${sample_id}_GATBMiniaPipeline --no-scaffolding --no-error-correction
        fi

        link=\$(readlink ${sample_id}_final.contigs.fa) && rm ${sample_id}_final.contigs.fa && mv \$link ${sample_id}_GATBMiniaPipeline.fasta

        echo pass > .status
    } || {
        echo fail > .status
        :> ${sample_id}_GATBMiniaPipeline.fasta
    }
    # rm temp dirs
    rm -r *_GATBMiniaPipeline.lib* *_GATBMiniaPipeline_besst *.unitigs* *contigs.fa *.h5 || true
    rm *list_reads* || true
    """
}


//      MINIA
IN_MINIA_kmer = Channel.value(params.miniakmer)

process MINIA {
    tag {sample_id}
    publishDir "results/$sample_id/assembly/MINIA/"

    input:
    tuple sample_id, file(fastq) from IN_MINIA
    val kmer from IN_MINIA_kmer

    output:
    tuple sample_id, val("MINIA"), file('*_minia.fasta') into OUT_MINIA
    file(".*version") into MINIA_VERSION

    script:
    """
    minia -v | head -n 1 | awk -F ' ' '{print \$3}' | awk -F 'v' '{print \$2}' | awk NF > .${sample_id}_MINIA_version
    {
        ls -1 $fastq  > list_reads
        minia -in list_reads -out ${sample_id}_minia.fasta -nb-cores $task.cpu
        mv ${sample_id}_minia.fasta.contigs.fa ${sample_id}_minia.fasta
        echo pass > .status
    } || {
        echo fail > .status
        :> ${sample_id}_minia.fasta
    }
    rm list_reads *.unitigs.* *.h5 || true
    """
}


//      MEGAHIT
IN_megahit_kmers = Channel.value(params.megahitKmers)

process MEGAHIT {
    tag { sample_id }
    publishDir "results/$sample_id/assembly/MEGAHIT/", pattern: '*_megahit*.fasta'

    input:
    tuple sample_id, file(fastq_pair) from IN_MEGAHIT
    val kmers from IN_megahit_kmers

    output:
    tuple sample_id, val("MEGAHIT"), file('*_MEGAHIT.fasta') into OUT_MEGAHIT
    file(".*version") into MEGAHIT_VERSION

    script:
    """
    /NGStools/megahit/bin/megahit -v | awk -F ' ' '{print \$2}' | awk -F 'v' '{print \$2}' | awk NF > .${sample_id}_MEGAHIT_version
    {
        /NGStools/megahit/bin/megahit --num-cpu-threads $task.cpus -o megahit --k-list $kmers -1 ${fastq_pair[0]} -2 ${fastq_pair[1]}
        mv megahit/final.contigs.fa ${sample_id}_MEGAHIT.fasta
        echo pass > .status
    } || {
        echo fail > .status
        :> ${sample_id}_MEGAHIT.fasta
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
    publishDir "results/$sample_id/assembly/metaSPAdes/"

    input:
    tuple sample_id, file(fastq_pair) from IN_METASPADES
    val kmers from IN_metaspades_kmers

    output:
    tuple sample_id, val("metaSPAdes"), file('*_metaspades.fasta') into OUT_METASPADES
    file(".*version") into METASPADES_VERSION

    script:
    """
    metaspades.py --version &> version
    cat version | awk -F ' ' '{print \$4}' | awk -F 'v' '{print \$2}' > .${sample_id}_metaSPAdes_version
    rm version
    {
        metaspades.py --only-assembler --threads $task.cpus -k $kmers -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} -o metaspades
        mv metaspades/contigs.fasta ${sample_id}_metaspades.fasta
        echo pass > .status
    } || {
        echo fail > .status
        :> ${sample_id}_metaspades.fasta
    }
    rm -r metaspades || true
    """
}


//      UNICYCLER
process UNICYCLER {
    tag { sample_id }
    publishDir "results/$sample_id/assembly/unicycler"

    input:
    tuple sample_id, file(fastq_pair) from IN_UNICYCLER

    output:
    tuple sample_id, val("Unicycler"), file('*_unicycler.fasta') into OUT_UNICYCLER
    file(".*version") into UNICYCLER_VERSION

    script:
    """
    unicycler --version | awk -F ' v' '{print \$2}' | awk NF > .${sample_id}_Unicycler_version 
    {
        unicycler -t $task.cpus -o . --no_correct --no_pilon -1 ${fastq_pair[0]} -2 ${fastq_pair[1]}
        mv assembly.fasta ${sample_id}_unicycler.fasta
        echo pass > .status
    } || {
        echo fail > .status
        :> ${sample_id}_unicycler.fasta
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
    publishDir "results/$sample_id/assembly/SPAdes/", pattern: '*.fasta'

    input:
    tuple sample_id, file(fastq_pair) from IN_SPADES
    val kmers from IN_spades_kmers

    output:
    tuple sample_id, val("SPAdes"), file('*_spades.fasta') into OUT_SPADES
    file(".*version") into SPADES_VERSION

    script:
    """
    spades.py --version &> version
    cat version | awk -F ' ' '{print \$4}' | awk -F 'v' '{print \$2}' > .${sample_id}_SPAdes_version
    rm version
    {
        spades.py --only-assembler --threads $task.cpus -k $kmers -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} -o spades
        mv spades/contigs.fasta ${sample_id}_spades.fasta
    } || {
        echo fail > .status
        :> ${sample_id}_spades.fasta
    }
    rm -r spades || true
    """
}


//      SKESA
process SKESA {
    tag { sample_id }
    publishDir "results/$sample_id/assembly/SKESA/"

    input:
    tuple sample_id, file(fastq_pair) from IN_SKESA

    output:
    tuple sample_id, val("SKESA"), file('*_skesa.fasta') into OUT_SKESA
    file(".*version") into SKESA_VERSION

    script:
    """
    skesa -v | tail -n 1 | awk -F ' ' '{print \$2}' | awk NF > .${sample_id}_SKESA_version
    {
        skesa --cores $task.cpus --memory $task.memory --use_paired_ends --contigs_out ${sample_id}_skesa.fasta --fastq ${fastq_pair[0]} ${fastq_pair[1]}
        echo pass > .status
    } || {
        echo fail > .status
        :> ${sample_id}_skesa.fasta
    }
    """
}


//      VELVETOPTIMIZER
process VELVETOPTIMIZER {
    tag { sample_id }
    publishDir "results/$sample_id/assembly/VelvetOtimiser"

    input:
    tuple sample_id, file(fastq_pair) from IN_VELVETOPTIMIZER

    output:
    tuple sample_id, val("VelvetOptimizer"), file('*.fasta') into OUT_VELVETOPTIMIZER
    file(".*version") into VELVETOPTIMIZER_VERSION

    script:
    """
    VelvetOptimiser.pl --version | awk -F ' ' '{print \$2}' | awk NF > .${sample_id}_VelvetOptimiser_version
    {
        VelvetOptimiser.pl -v -s $params.velvetoptimizer_hashs -e $params.velvetoptimizer_hashe -t $task.cpus \
        -f '-shortPaired -fastq.gz -separate ${fastq_pair[0]} ${fastq_pair[1]}'
        mv auto_data*/contigs.fa ${sample_id}_velvetoptimizer.fasta
        echo pass > .status
    } || {
        echo fail > .status
        :> ${sample_id}_velvetoptimizer.fasta
    }
    rm -r auto_data* || true
    """
}


//      IDBA
process reformat_IDBA {
    tag { sample_id }

    input:
    tuple sample_id, file(fastq_pair) from IN_IDBA

    output:
    tuple sample_id, file('*.fasta') into REFORMAT_IDBA

    script:
    "reformat.sh in=${fastq_pair[0]} in2=${fastq_pair[1]} out=${sample_id}_reads.fasta"
}

process IDBA {
    tag { sample_id }
    publishDir "results/$sample_id/assembly/IDBA-UD/"

    input:
    tuple sample_id, file(fasta_reads_single) from  REFORMAT_IDBA

    output:
    tuple sample_id, val("IDBA-UD"), file('*_IDBA-UD.fasta') into OUT_IDBA
    file(".*version") into IDBA_VERSION

    script:
    """
    echo '' > .${sample_id}_IDBA_version
    {
        idba_ud -l ${fasta_reads_single} --num_threads $task.cpus -o .
        mv contig.fa ${sample_id}_IDBA-UD.fasta
        echo pass > .status
    } || {
        echo fail > .status
        :> ${sample_id}_IDBA-UD.fasta
    }
    rm begin align-* contig-* graph-* kmer local-* || true
    """
}

// VERSION COLLECTION
BCALM2_VERSION.mix(GATB_VERSION,
                    MINIA_VERSION,
                    MEGAHIT_VERSION,
                    METASPADES_VERSION,
                    UNICYCLER_VERSION,
                    SPADES_VERSION,
                    SKESA_VERSION,
                    VELVETOPTIMIZER_VERSION,
                    IDBA_VERSION).set{ALL_VERSIONS}

process PROCESS_VERSION {

    input:
    file version from ALL_VERSIONS.collect()

    output:
    file("versions.json") into VERSIONS_JSON

    script:
    template "process_versions.py"
}

// ASSEMBLY COLLECTION
OUT_BCALM2.mix(OUT_GATB,
                  OUT_MINIA,
                  OUT_MEGAHIT,
                  OUT_METASPADES,
                  OUT_UNICYCLER,
                  OUT_SPADES,
                  OUT_SKESA,
                  OUT_VELVETOPTIMIZER,
                  OUT_IDBA).set{ALL_ASSEMBLERS}

ALL_ASSEMBLERS.into{ TO_FILTER; TO_GLOBAL_STATS; TO_READ_MAPPING}

THRESHOLD = Channel.value(params.mapped_reads_threshold)
// READ MAPPING
process READ_MAPPING{

    tag { assembler }

    publishDir "results/$sample_id/mapping/reads"

    input:
    tuple sample_id, assembler, assembly from TO_READ_MAPPING
    each THRESHOLD

    output:
    file("*_read_mapping.txt") optional true
    tuple sample_id, assembler, file("*_read_mapping_report.json") into OUT_READ_MAPPING optional true

    script:
    template "read_mapping.py"
}

// ASSEMBLY STATS GLOBAL
process ASSEMBLY_STATS_GLOBAL {
    tag { assembler }

    publishDir "results/$sample_id/stats/assembly"

    input:
    tuple sample_id, assembler, file(assembly), file(read_mapping) from TO_GLOBAL_STATS.join(OUT_READ_MAPPING, by: [0,1])

    output:
    file "*report.json" into OUT_ASSEMBLY_STATS_GLOBAL_JSON
    file "*.csv" into OUT_ASSEMBLY_STATS_GLOBAL_TSV

    script:
    template "assembly_stats_global.py"
}


process PROCESS_ASSEMBLY_STATS_GLOBAL {

    publishDir "results/stats"

    input:
    file assembly_stats_global_files from OUT_ASSEMBLY_STATS_GLOBAL_TSV.collect()
    file json_report from OUT_ASSEMBLY_STATS_GLOBAL_JSON.collect()

    output:
    file "global_assembly_stats.json" into PROCESS_ASSEMBLY_STATS_GLOBAL_OUT

    script:
    template "process_assembly_stats_global.py"

}

// FILTER ASSEMBLY
IN_minLen = Channel.value(params.minLength)

process FILTER_ASSEMBLY {

    tag {sample_id; assembler}
    publishDir "results/$sample_id/assembly/filtered/"

    input:
    tuple sample_id, assembler, file(assembly) from TO_FILTER
    val minLen from IN_minLen

    output:
    tuple sample_id, assembler, file('filtered_*') into OUT_FILTERED

    script:
    "reformat.sh in=${assembly} out=filtered_${assembly} minlength=${minLen}"
}

// ASSEMBLY MAPPING
process ASSEMBLY_MAPPING{

    tag { sample_id; assembler }

    publishDir "results/$sample_id/mapping/assembly"

    input:
    tuple sample_id, assembler, file(assembly) from OUT_FILTERED
    each reference from IN_MAPPING_CONTIGS

    output:
    tuple sample_id, assembler, file(assembly), file("*.paf") into OUT_ASSEMBLY_MAPPING

    script:
    "minimap2 --cs -N 0 -t $task.cpus -r 10000 -g 10000 -x asm20 --eqx ${reference} ${assembly} > ${sample_id}_${assembler}.paf"

}

OUT_ASSEMBLY_MAPPING.into{ IN_ASSEMBLY_MAPPING_FOR_STATS; IN_GAP_ASSESSMENT; IN_SNP_ASSESSMENT; IN_MISASSEMBLY}

process ASSEMBLY_STATS_MAPPING {

    tag { assembler }

    publishDir "results/$sample_id/stats/"

    input:
    tuple sample_id, assembler, file(assembly), file(mapping) from IN_ASSEMBLY_MAPPING_FOR_STATS
    each reference from IN_ASSEMBLY_STATS_MAPPING

    output:
    file("*_report.json") into OUT_ASSEMBLY_STATS_MAPPING_JSON
    file("*breadth_of_coverage_contigs.csv") into OUT_COVERAGE_PER_CONTIG
    file "*_df.csv" into OUT_DF_ASSEMBLY_STATS_MAPPING
    file("*_lx.csv") into OUT_LX_PLOT
    file("*_nax.csv") into OUT_NAX_PLOT
    file("*_ngx.csv") into OUT_NGX_PLOT
    file("*_phred.csv") into OUT_PHRED

    script:
    template "assembly_stats_mapping.py"

}

process PROCESS_ASSEMBLY_STATS_MAPPING {

    publishDir "results/stats/"

    input:
    file json_report from OUT_ASSEMBLY_STATS_MAPPING_JSON.collect()

    output:
    file "global_assembly_mapping_stats.json" into PROCESS_ASSEMBLY_STATS_MAPPING_OUT

    script:
    template "process_assembly_stats_mapping.py"

}

process PROCESS_COMPLETNESS {

    publishDir 'results/plots/', pattern: "*.html", mode: "copy"

    input:
    file coverage_files from OUT_COVERAGE_PER_CONTIG.collect()

    output:
    file("*.html")
    file("completness_plots.json") into PLOT_PROCESS_COMPLETNESS

    script:
    template "completness_plot.py"
}

process PLOT_LX {

    publishDir 'results/plots/', pattern: "*.html", mode: "copy"

    input:
    file lx_files from OUT_LX_PLOT.collect()

    output:
    file("*.html")
    file("lx.json") into PLOT_LX

    script:
    template "lx_plot.py"
}

process PLOT_NAX {

    publishDir 'results/plots/', pattern: "*.html", mode: "copy"

    input:
    file nax_files from OUT_NAX_PLOT.collect()

    output:
    file("*.html")
    file("nax.json") into PLOT_NAX

    script:
    template "nax_plot.py"
}

process PLOT_NGX {

    publishDir 'results/plots/', pattern: "*.html", mode: "copy"

    input:
    file ngx_files from OUT_NGX_PLOT.collect()

    output:
    file("*.html")
    file("ngx.json") into PLOT_NGX

    script:
    template "ngx_plot.py"
}

process PROCESS_SHRIMP_PLOT {

    publishDir 'results/plots/', pattern: "*.html", mode: "copy"

    input:
    file phred_files from OUT_PHRED.collect()

    output:
    file("*.html")
    file("phred.json") into PLOT_PHRED

    script:
    template "shrimp_plot.py"
}


process PLOT_CONTIG_DISTRIBUTION {

    publishDir 'results/plots/', pattern: "*.html", mode: "copy"

    input:
    file dataframes from OUT_DF_ASSEMBLY_STATS_MAPPING.collect()

    output:
    file("*.html")
    file("*.json") into PLOT_CONTIG_DISTRIBUTION

    script:
    template "plot_contig_size.py"
}

process GAP_ASSESSMENT {

    tag { assembler }

    input:
    tuple sample_id, assembler, file(assembly), file(mapping) from IN_GAP_ASSESSMENT
    each reference from IN_GAP_STATS

    output:
    file("*_gap_dict.json") into OUT_GAP_DISTANCE
    file("*_gaps.csv") into OUT_GAP_PLOT_REF

    script:
    template "gap_assessment.py"
}

process PLOT_GAP_BOXPLOT {

    publishDir 'results/plots/', pattern: "*.html", mode: "copy"

    input:
    file gap_distance_json from OUT_GAP_DISTANCE.collect()

    output:
    file("*.html")
    file("*gap_distance_histogram.json") into OUT_GAP_HISTOGRAM

    script:
    template "plot_gap_sizes.py"

}

process PLOT_GAP_REFERENCE {

    publishDir 'results/plots/', pattern: "*.html", mode: "copy"

    input:
    file gap_coords_dataframes from OUT_GAP_PLOT_REF.collect()

    output:
    file("*.html")
    file("*.json") into OUT_GAP_REFERENCE

    script:
    template "plot_gap_reference.py"
}

process SNP_ASSESSMENT {

    tag { assembler }

    input:
    tuple sample_id, assembler, file(assembly), file(mapping) from IN_SNP_ASSESSMENT
    each reference from IN_SNP_STATS

    output:
    file("*.tsv")
    file("*_snps.csv") into OUT_SNP_PLOT_REF

    script:
    template "snp_assessment.py"
}


process PLOT_SNP_REFERENCE {

    publishDir 'results/plots/', pattern: "*.html", mode: "copy"

    input:
    file snp_coords_dataframes from OUT_SNP_PLOT_REF.collect()

    output:
    file("*.html")
    file("*.json") into OUT_SNP_REFERENCE

    script:
    template "plot_snp.py"
}

process MISASSEMBLY {

    tag { assembler }

    input:
    tuple sample_id, assembler, file(assembly), file(mapping) from IN_MISASSEMBLY

    output:
    file("*_trace.pkl") into OUT_MISASSEMBLY_TRACE
    file("*_contig_lenght.pkl") into OUT_MISASSEMBLY_CONTIGS
    file("*_misassembly.json") into MISASSEMBLY_REPORT
    file("*_misassembled_reference.json") into MISASSEMBLY_DICTIONARY
    file("*_misassembly.csv") into PLOT_MISASSEMBLY_REF

    script:
    template "misassembly.py"

}

process PROCESS_MISASSEMBLY {

    publishDir 'results/plots/', pattern: "*.html", mode: "copy"

    input:
    file misassembly_trace from OUT_MISASSEMBLY_TRACE.collect()
    file misassembly_contigs from OUT_MISASSEMBLY_CONTIGS.collect()
    file report_data from MISASSEMBLY_REPORT.collect()
    file report_per_reference from MISASSEMBLY_DICTIONARY.collect()

    output:
    file("*.html")
    file("*_misassembly.json") into OUT_MISASSEMBLY_PLOT
    file("misassembly_report.json") into OUT_MISASSEMBLY_REPORT
    file("misassembly_report_per_ref.json") into MISASSEMBLY_PER_REF

    script:
    template "process_misassembly.py"

}

process PLOT_MISASSEMBLY {

    publishDir 'results/plots/', pattern: "*.html", mode: "copy"

    input:
    file misassembly_dataframes from PLOT_MISASSEMBLY_REF.collect()

    output:
    file("*.html")
    file("*.json") into OUT_MISASSEMBLY_REFERENCE

    script:
    template "plot_misassembly.py"


}

/** Reports
Compiles the reports from every process
**/

OUT_ASSEMBLY_STATS_GLOBAL_JSON.set{master_report}

process compile_reports {

    publishDir "report/", mode: "copy"

    input:
    file reads_json from PROCESS_READS.collect()
    file global_assembly_stats from PROCESS_ASSEMBLY_STATS_GLOBAL_OUT
    file pipeline_stats from Channel.fromPath("${workflow.projectDir}/pipeline_stats.txt")
    file js from Channel.fromPath("${workflow.projectDir}/resources/main.js.zip")
    file lmas_png from Channel.fromPath("${workflow.projectDir}/resources/lmas.zip")
    file reference_file from COMPILE_REPORTS_REF
    file contig_size_distribution from PLOT_CONTIG_DISTRIBUTION
    file mapping_assembly_stats from PROCESS_ASSEMBLY_STATS_MAPPING_OUT
    file completness_plots from PLOT_PROCESS_COMPLETNESS
    file lx_plots from PLOT_LX
    file shrimp_plots from PLOT_PHRED
    file gap_reference_json from OUT_GAP_REFERENCE
    file snp_reference_json from OUT_SNP_REFERENCE
    file gap_histogram from OUT_GAP_HISTOGRAM
    file plot_misassemblies from OUT_MISASSEMBLY_PLOT
    file misassembly_data from OUT_MISASSEMBLY_REPORT
    file nax_plots from PLOT_NAX
    file ngx_plots from PLOT_NGX
    file versions_json from VERSIONS_JSON
    file misassembly_per_ref from MISASSEMBLY_PER_REF
    file plot_misassembly_per_ref from OUT_MISASSEMBLY_REFERENCE
    file about_md from Channel.fromPath(params.md)
    file containers_config from Channel.fromPath("${workflow.projectDir}/containers.config")

    output:
    file "pipeline_report*.json"
    file "index.html"
    file "main.js"
    file "lmas.jpg"
    file "performance_metadata.json"
    file "reference_metadata.json"

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