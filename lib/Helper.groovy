class Help {

    static def start_info(Map info, String time, String profile, String version) {

        println ""
        println "                      _    __  __   _   ___"
        println "       /\\︵︵/\\      | |  |  \\/  | /_\\ / __|"
        println "      (◕('人')◕)     | |__| |\\/| |/ _ \\\\__ \\"
        println "         |︶|        |____|_|  |_/_/ \\_\\___/"
        println ""
        println "         Last Metagenomic Assembler Standing"
        println ""
        int nsamples = info.fastq / 2
        println " Input FastQ                 : $nsamples"
        println " Reference file              : $info.referece"
        println " Reports are found in        : ./reports"
        println " Results are found in        : ./results"
        println " Profile                     : $profile"
        println " Version                     : $version"
        println ""
        println "Starting pipeline at $time"
        println ""

    }

    static void complete_info(nextflow.script.WorkflowMetadata wf) {

        println ""
        println "Pipeline execution summary"
        println "=========================="
        println "Completed at                 : $wf.complete"
        println "Duration                     : $wf.duration"
        println "Success                      : $wf.success"
        println "Work directory               : $wf.workDir"
        println "Exit status                  : $wf.exitStatus"
        println ""

    }

    static def print_help(Map params) {

        println ""
        println "                      _    __  __   _   ___"
        println "       /\\︵︵/\\      | |  |  \\/  | /_\\ / __|"
        println "      (◕('人')◕)     | |__| |\\/| |/ _ \\\\__ \\"
        println "         |︶|        |____|_|  |_/_/ \\_\\___/"
        println ""
        println "         Last Metagenomic Assembler Standing"
        println ""
        println "Basic Usage: "
        println "    nextflow run LMAS.nf"
        println ""
        println "Input parameters:"
        println "    --fastq                    Path expression to paired-end fastq files." 
        println "                               (default: $params.fastq)"
        println "    --reference                Path to the genome reference fasta file."
        println "                               (default: $params.reference)"
        println "    --md                       Path to markdown with input sample description for report (optional)."
        println "                               (default: $params.md)"
        println ""
        println "Mapping and filtering paramenters:"
        println "    --minLength                Value for minimum contig length, in basepairs."
        println "                               (default: $params.minLength)"
        println "    --mapped_reads_threshold   Value for the minimum percentage of a read aligning to the"
        println "                               contig to be considered as mapped."
        println "                               (default: $params.mapped_reads_threshold)"
        println ""
        println "Assembly quality assessment parameters:"
        println "    --n_target                 Target value for the N, NA and NG metrics, ranging from 0 to 1."
        println "                               (default: $params.n_target)"
        println "    --l_target                 Target value for the L metric, ranging from 0 to 1."
        println "                               (default: $params.n_target)"
        println "    --plot_scale               Scale of x-axis for the L, NA and NG metrics plots."
        println "                               Allowed values: 'linear' or 'log'."
        println "                               (default: $params.plot_scale)"
        println ""
        println "Assembly execution parameters:"
        println "    --abyss                    Boolean controling the execution of the ABySS assembler."
        println "                               (default: $params.abyss)"
        println "    --abyssKmerSize            K-mer size for the ABySS assembler, as an intiger."
        println "                               (default $params.abyssKmerSize)"
        println "    --abyssBloomSize           Bloom filter size for the ABySS assembler."
        println "                               It must be a sting with a value and an unit."
        println "                               (default: $params.abyssBloomSize)"
        println "    --bcalm                    Boolean controling the execution of the BCALM2 assembler."
        println "                               (default: $params.bcalm)"
        println "    --bcalmKmerSize            K-mer size for the BCALM2 assembler, as an intiger."
        println "                               (default $params.bcalmKmerSize)"
        println "    --gatb_minia               Boolean controling the execution of the GATB Minia Pipeline assembler."
        println "                               (default: $params.gatb_minia)"
        println "    --gatbKmerSize             K-mer sizes for the GATB Minia Pipeline assembler."
        println "                               It must be a sting with the values separated with a comma."
        println "                               (default $params.gatbKmerSize)"
        println "    --gatb_besst_iter          Number of iteration during Besst scaffolding for the"
        println "                               GATB Minia Pipeline assembler."
        println "                               (default $params.gatb_besst_iter)"
        println "    --gatb_error_correction    Boolean to control weather to skip error correction for the"
        println "                               GATB Minia Pipeline assembler."
        println "                               (default $params.gatb_error_correction)"
        println "    --idba                     Boolean controling the execution of the IDBA-UD assembler."
        println "                               (default $params.idba)"
        println "    --metahipmer2              Boolean controling the execution of the MetaHipMer2 assembler."
        println "                               (default $params.metahipmer2)"
        println "    --metahipmer2KmerSize      K-mer sizes for the MetaHipMer2 assembler."
        println "                               It must be a sting with the values separated with a comma."
        println "                               (default $params.metahipmer2KmerSize)"
        println "    --minia                    Boolean controling the execution of the minia assembler."
        println "                               (default: $params.minia)"
        println "    --miniaKmerSize            K-mer size for the minia assembler, as an intiger."
        println "                               (default $params.miniaKmerSize)"
        println "    --megahit                  Boolean controling the execution of the MEGAHIT assembler."
        println "                               (default $params.megahit)"
        println "    --megahitKmerSize          K-mer sizes for the MEGAHIT assembler."
        println "                               It must be a sting with the values separated with a comma."
        println "                               (default $params.megahitKmerSize)"
        println "    --metaspades               Boolean controling the execution of the metaSPAdes assembler."
        println "                               (default $params.metaspades)"
        println "    --metaspadesKmerSize       K-mer sizes for the metaSPAdes assembler."
        println "                               It must be a sting with 'auto' or the values separated with a space."
        println "                               (default $params.metaspadesKmerSize)"
        println "    --spades                   Boolean controling the execution of the SPAdes assembler."
        println "                               (default $params.spades)"
        println "    --spadesKmerSize           K-mer sizes for the SPAdes assembler."
        println "                               It must be a sting with 'auto' or the values separated with a space."
        println "                               (default $params.spadesKmerSize)"
        println "    --skesa                    Boolean controling the execution of the SKESA assembler."
        println "                               (default $params.skesa)"
        println "    --unicycler                Boolean controling the execution of the Unicycler assembler."
        println "                               (default $params.unicycler)"
        println "    --velvetoptimiser          Boolean controling the execution of the VelvetOptimiser assembler."
        println "                               (default: $params.velvetoptimiser)"
        println "    --velvetoptimiser_hashs    Starting K-mer size for the VelvetOptimiser assembler, as an intiger."
        println "                               (default $params.velvetoptimiser_hashs)"
        println "    --velvetoptimiser_hashe    End K-mer size for the VelvetOptimiser assembler, as an intiger."
        println "                               (default $params.velvetoptimiser_hashe)"
        println ""
        println "Execution resources parameters:"
        println "    --cpus                     Number of CPUs for the assembly and mapping processes, as an intiger."
        println "                               This resource is double for each retry until max_cpus is reached."
        println "                               (default $params.cpus)"
        println "    --memory                   Memory for the assembly and mapping processes, in the format of" 
        println "                               'value'.'unit'."
        println "                               This resource is double for each retry until max_memory is reached."
        println "                               (default $params.memory)"
        println "    --time                     Time limit for the assembly and mapping processes, in the format of" 
        println "                               'value'.'unit'."
        println "                               This resource is double for each retry until max_time is reached."
        println "                               (default $params.time)"
        println "    --max_cpus                 Maximum number of CPUs for the assembly and mapping processes,"
        println "                               as an intiger. It overwrites the --cpu parameter."
        println "                               (default $params.max_cpus)"
        println "    --max_memory               Maximum memory for the assembly and mapping processes, in the format of"
        println "                               'value'.'unit'. It overwrites the --memory parameter."
        println "                               (default $params.max_memory)"
        println "    --max_time                 Maximum time for the assembly and mapping processes, in the format of"
        println "                               'value'.'unit'. It overwrites the --time parameter."
        println "                               (default $params.max_time)"
    }

}

class CollectInitialMetadata {

    public static void print_metadata(nextflow.script.WorkflowMetadata workflow){

        def metadataJson = "{'nfMetadata':{'scriptId':'${workflow.scriptId}',\
                            'scriptName':'${workflow.scriptName}',\
                            'profile':'${workflow.profile}',\
                            'container':'${workflow.container}',\
                            'containerEngine':'${workflow.containerEngine}',\
                            'commandLine':'${workflow.commandLine}',\
                            'runName':'${workflow.runName}',\
                            'sessionId':'${workflow.sessionId}',\
                            'projectDir':'${workflow.projectDir}',\
                            'launchDir':'${workflow.launchDir}',\
                            'startTime':'${workflow.start}'}}"

        def json = metadataJson.replaceAll("'", '"')

        def jsonFile = new File(".metadata.json")
        jsonFile.write json
    }
}
