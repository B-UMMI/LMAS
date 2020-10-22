class Help {

    static def start_info(Map info, String time, String profile, String version) {

        println ""
        println "LMAS - Last Metagenomic Assembler Standing"
        println ""
        int nsamples = info.fastq / 2
        println " Input FastQ                 : $info.fastq"
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
        println "LMAS - Last Metagenomic Assembler Standing"
        println ""
        println "Usage: "
        println "    nextflow run LMAS.nf"
        println ""
        println "    --fastq               Path expression to fastq files. (default: $params.fastq)"
        println ""
        println"     --reference           Path to triple-genome reference fasta file. (default: $params.referece)"
        println ""

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
