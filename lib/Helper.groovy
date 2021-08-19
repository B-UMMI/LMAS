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
        println "    --reference                Path to triple-genome reference fasta file."
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
