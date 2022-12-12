class Params {

    static void check(Map params) {

        // mode
        if (!params.wf) {print_error("'--wf' parameter missing")}
        def allowed_modes = ["default", "Illumina", "illumina", "hybrid", "Hybrid", "ont", "ONT"] as Set
        def mode_parameter_diff = allowed_modes - params.wf
        if (mode_parameter_diff.size() > 6){
                print_error("[Pipeline warning] Parameter --wf $params.wf is not valid in the pipeline! Allowed values: 'default', 'long', 'hybrid', 'all'")
    }
        
        // input
        if (!params.reference) {print_error("'--reference' parameter missing")}
        if (!params.fastq) {print_error("'--fastq' parameter missing")}
        if (params.reference instanceof Boolean) {print_error("'--reference' must be a path pattern. Provided value: '$params.reference'")}
        if (params.fastq instanceof Boolean) {print_error("'--fastq' must be a path pattern. Provided value: '$params.fastq'")}

        // assembler skipping
        if (!params.abyss && !params.gatb_minia && !params.idba && !params.metahipmer2 && !params.minia && !params.megahit && !params.metaspades && !params.spades && !params.skesa && !params.unicycler && !params.velvetoptimiser) {print_error("All assembly processes set to false. Exiting.")}

        // assembler parameters
        if (!params.abyssKmerSize.toString().isNumber()) {print_error("'--abyssKmerSize' parameter must be a number. Provided value: '$params.abyssKmerSize'")}
        if (!params.gatb_besst_iter.toString().isNumber()) {print_error("'--gatb_besst_iter' parameter must be a number. Provided value: '$params.gatb_besst_iter'")}
        if (params.metaspadesKmerSize.toString().split(" ").size() <= 1) {if (params.metaspadesKmerSize.toString() != 'auto') {print_error("'--metaspadesKmerSize' parameter must be a sequence of space separated numbers or 'auto'. Provided value: '$params.metaspadesKmerSize'")}}
        if (params.spadesKmerSize.toString().split(" ").size() <= 1) {if (params.spadesKmerSize.toString() != 'auto'){print_error("'--spadesKmerSize' parameter must be a sequence of space separated numbers or 'auto'. Provided value: '$params.spadesKmerSize'")}}
        if (!params.minLength.toString().isNumber()) {print_error("'--minLength' parameter must be a number. Provided value: '$params.minLength'")}

    }

    static def print_error(String msg) {

        println "\nERROR: $msg"
        System.exit(1)

    }

}
