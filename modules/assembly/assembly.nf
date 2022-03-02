nextflow.enable.dsl=2

// PROCESSES
process REFORMAT {

    tag { sample_id }
    label 'process_assembly'

    when:
    params.idba || params.metahipmer2

    input:
    tuple val(sample_id), path(fastq_pair)

    output:
    tuple val(sample_id), file('*.fasta') 

    script:
    "reformat.sh in=${fastq_pair[0]} in2=${fastq_pair[1]} out=${sample_id}_reads.fasta"
}

process ABYSS {

    tag { sample_id }
    label 'process_assembly'
    publishDir "results/$sample_id/assembly/abyss/"

    when:
    params.abyss

    input:
    tuple val(sample_id), path(fastq)
    val KmerSize 
    val BloomSize 

    output:
    tuple val(sample_id), val('ABySS'), path('*_ABySS.fasta'), emit: assembly
    path('.*version'), emit: version

    script:
    """
    abyss-pe version | grep "ABySS" | awk -F ' ' '{print \$3}' > .${sample_id}_ABySS_version
    {
        abyss-pe name='${sample_id}' j=$task.cpus k=$KmerSize B=$BloomSize in='$fastq'
        mv ${sample_id}-contigs.fa ${sample_id}_ABySS.fasta
        echo pass > .status
    } || {
        echo fail > .status
        :> ${sample_id}_ABySS.fasta
    }
    # remove temp files
    rm  *.dot* *.hist *.path* || true
    """
}

process BCALM2 {

    tag { sample_id }
    label 'process_assembly'
    publishDir "results/$sample_id/assembly/bcalm2/"

    when:
    params.bcalm

    input:
    tuple val(sample_id), path(fastq) 
    val KmerSize 

    output:
    tuple val(sample_id), val('BCALM2'), path('*_BCALM2.fasta'), emit: assembly
    path('.*version'), emit: version 

    script:
    """
    ls -1 $fastq  > list_reads
    bcalm -version | head -n 1 | awk -F ', ' '{print \$2}' | awk -F ' ' '{print \$2}' | awk -F 'v' '{print \$2}' \
    > .${sample_id}_BCALM2_version
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

process GATBMINIAPIPELINE {

    tag { sample_id }
    label 'process_assembly'
    publishDir "results/$sample_id/assembly/GATBMiniaPipeline/"

    when:
    params.gatb_minia

    input:
    tuple val(sample_id), path(fastq_pair)
    val kmer_list
    val do_error_correction
    val besst_iter

    output:
    tuple val(sample_id), val('GATBMiniaPipeline'), path('*_GATBMiniaPipeline.fasta'), emit: assembly
    path('.*version'), emit: version 

    script:
    """
    echo '' > .${sample_id}_GATBMiniaPipeline_version
    {
        if [ $do_error_correction ];
        then
            gatb -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} --kmer-sizes ${kmer_list} \
            -o ${sample_id}_GATBMiniaPipeline --no-scaffolding
        else
            gatb -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} --kmer-sizes ${kmer_list} \
            -o ${sample_id}_GATBMiniaPipeline --no-scaffolding --no-error-correction
        fi

        link=\$(readlink *_final.contigs.fa) && mv \$link ${sample_id}_GATBMiniaPipeline.fasta

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

process IDBA {

    tag { sample_id }
    label 'process_assembly'
    publishDir "results/$sample_id/assembly/IDBA-UD/"

    when:
    params.idba

    input:
    tuple val(sample_id), path(fasta_reads_single)

    output:
    tuple val(sample_id), val('IDBA-UD'), file('*_IDBA-UD.fasta'), emit: assembly
    path('.*version'), emit: version 

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

process MEGAHIT {

    tag { sample_id }
    label 'process_assembly'
    publishDir "results/$sample_id/assembly/MEGAHIT/", pattern: '*.fasta'

    when:
    params.megahit

    input:
    tuple val(sample_id), path(fastq_pair)
    val kmers

    output:
    tuple val(sample_id), val('MEGAHIT'), path('*_MEGAHIT.fasta'), emit: assembly
    path('.*version'), emit: version

    script:
    """
    /NGStools/megahit/bin/megahit -v | awk -F ' ' '{print \$2}' | awk -F 'v' '{print \$2}' | awk NF \
    > .${sample_id}_MEGAHIT_version
    {
        /NGStools/megahit/bin/megahit --num-cpu-threads $task.cpus -o megahit --k-list $kmers \
        -1 ${fastq_pair[0]} -2 ${fastq_pair[1]}

        mv megahit/final.contigs.fa ${sample_id}_MEGAHIT.fasta
        echo pass > .status
    } || {
        echo fail > .status
        :> ${sample_id}_MEGAHIT.fasta
    }
    rm -r megahit || true
    """
}

process METAHIPMER2 {

    tag { sample_id }
    label 'process_assembly'
    publishDir "results/$sample_id/assembly/MetaHipMer2/"

    when:
    params.metahipmer2

    input:
    tuple val(sample_id), path(fasta_reads_single)
    val kmer

    output:
    tuple val(sample_id), val('MetaHipMer2'), file('*_MetaHipMer2.fasta'), emit: assembly
    path('.*version'), emit: version

    script:
    """
    mhm2.py -h | grep "version" |  awk -F ' ' '{print \$3}' > .${sample_id}_MetaHipMer2_version
    {
        mhm2.py -r $fasta_reads_single -k $kmer -s 0 --max-kmer-store 20 --procs $task.cpus  \
        --max-rpcs-in-flight 50 --shared-heap 800mb

        mv mhm2-run*/final_assembly.fasta ${sample_id}_MetaHipMer2.fasta
        echo pass > .status
    } || {
        echo fail > .status
        :> ${sample_id}_MetaHipMer2.fasta
    }
    rm -r mhm2-run* || true
    """
}

process METASPADES {

    tag { sample_id }
    label 'process_assembly'
    publishDir "results/$sample_id/assembly/metaSPAdes/"

    when:
    params.metaspades

    input:
    tuple val(sample_id), path(fastq_pair)
    val kmers

    output:
    tuple val(sample_id), val('metaSPAdes'), path('*_metaspades.fasta'), emit: assembly
    path('.*version'), emit: version

    script:
    """
    metaspades.py --version &> version
    cat version | awk -F ' ' '{print \$4}' | awk -F 'v' '{print \$2}' > .${sample_id}_metaSPAdes_version
    rm version
    {
        metaspades.py --only-assembler --threads $task.cpus -k $kmers \
        -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} -o metaspades

        mv metaspades/contigs.fasta ${sample_id}_metaspades.fasta
        echo pass > .status
    } || {
        echo fail > .status
        :> ${sample_id}_metaspades.fasta
    }
    rm -r metaspades || true
    """
}

process MINIA {

    tag {sample_id}
    label 'process_assembly'
    publishDir "results/$sample_id/assembly/MINIA/"

    when:
    params.minia

    input:
    tuple val(sample_id), path(fastq)
    val kmer

    output:
    tuple val(sample_id), val('MINIA'), path('*_minia.fasta'), emit: assembly
    path('.*version'), emit: version

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

process SKESA {

    tag { sample_id }
    label 'process_assembly'
    publishDir "results/$sample_id/assembly/SKESA/"

    when:
    params.skesa

    input:
    tuple val(sample_id), path(fastq_pair)

    output:
    tuple val(sample_id), val('SKESA'), path('*_skesa.fasta'), emit: assembly
    path('.*version'), emit: version 

    script:
    """
    skesa -v | tail -n 1 | awk -F ' ' '{print \$2}' | awk NF > .${sample_id}_SKESA_version
    {
        skesa --cores $task.cpus --memory $task.memory --use_paired_ends --contigs_out ${sample_id}_skesa.fasta \
        --fastq ${fastq_pair[0]} ${fastq_pair[1]}

        echo pass > .status
    } || {
        echo fail > .status
        :> ${sample_id}_skesa.fasta
    }
    """
}

process SPADES {

    tag { sample_id }
    label 'process_assembly'
    publishDir "results/$sample_id/assembly/SPAdes/", pattern: '*.fasta'

    when:
    params.spades

    input:
    tuple val(sample_id), file(fastq_pair) 
    val kmers

    output:
    tuple val(sample_id), val('SPAdes'), path('*_spades.fasta'), emit: assembly
    path('.*version'), emit: version

    script:
    """
    spades.py --version &> version
    cat version | awk -F ' ' '{print \$4}' | awk -F 'v' '{print \$2}' > .${sample_id}_SPAdes_version
    rm version
    {
        spades.py --only-assembler --threads $task.cpus -k $kmers \
        -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} -o spades

        mv spades/contigs.fasta ${sample_id}_spades.fasta
    } || {
        echo fail > .status
        :> ${sample_id}_spades.fasta
    }
    rm -r spades || true
    """
}

process UNICYCLER {

    tag { sample_id }
    label 'process_assembly'
    publishDir "results/$sample_id/assembly/unicycler"

    when:
    params.unicycler

    input:
    tuple val(sample_id), file(fastq_pair)

    output:
    tuple val(sample_id), val('Unicycler'), path('*_unicycler.fasta'), emit: assembly
    path('.*version'), emit: version

    script:
    """
    unicycler --version | awk -F ' v' '{print \$2}' | awk NF > .${sample_id}_Unicycler_version
    {
        unicycler -t $task.cpus -o . --no_correct --no_pilon \
        -1 ${fastq_pair[0]} -2 ${fastq_pair[1]}

        mv assembly.fasta ${sample_id}_unicycler.fasta
        echo pass > .status
    } || {
        echo fail > .status
        :> ${sample_id}_unicycler.fasta
    }
    rm *best_spades_graph* *overlaps_removed* *bridges_applied* *final_clean* || true
    """
}

process VELVETOPTIMISER {

    tag { sample_id }
    label 'process_assembly'
    publishDir "results/$sample_id/assembly/VelvetOptimiser"

    when:
    params.velvetoptimiser

    input:
    tuple val(sample_id), path(fastq_pair)

    output:
    tuple val(sample_id), val('VelvetOptimiser'), path('*.fasta'), emit: assembly
    path('.*version'), emit: version

    script:
    """
    VelvetOptimiser.pl --version | awk -F ' ' '{print \$2}' | awk NF > .${sample_id}_VelvetOptimiser_version
    {
        VelvetOptimiser.pl -v -s $params.velvetoptimiser_hashs -e $params.velvetoptimiser_hashe -t $task.cpus \
        -f '-shortPaired -fastq.gz -separate ${fastq_pair[0]} ${fastq_pair[1]}'

        mv auto_data*/contigs.fa ${sample_id}_velvetoptimiser.fasta
        echo pass > .status
    } || {
        echo fail > .status
        :> ${sample_id}_velvetoptimiser.fasta
    }
    rm -r auto_data* || true
    """
}

// WORKFLOWS
workflow assembly_wf {

    abyssKmerSize = Channel.value(params.abyssKmerSize)
    abyssBloomSize = Channel.value(params.abyssBloomSize)
    bcalmKmerSize = Channel.value(params.bcalmKmerSize)
    gatbKmerSize = Channel.value(params.gatbKmerSize)
    GATB_error_correction = params.gatb_error_correction ? 'true' : 'false'
    gatb_besst_iter = Channel.value(params.gatb_besst_iter)
    megahitKmerSize = Channel.value(params.megahitKmerSize)
    metahipmer2KmerSize = Channel.value(params.metahipmer2KmerSize)
    metaspadesKmerSize = Channel.value(params.metaspadesKmerSize)
    miniaKmerSize = Channel.value(params.miniaKmerSize)
    spadesKmerSize = Channel.value(params.spadesKmerSize)

    take:
    IN_fastq_raw

    main:
    REFORMAT(IN_fastq_raw)
    ABYSS(IN_fastq_raw, abyssKmerSize, abyssBloomSize)
    BCALM2(IN_fastq_raw, bcalmKmerSize)
    GATBMINIAPIPELINE(IN_fastq_raw, gatbKmerSize, GATB_error_correction, gatb_besst_iter)
    IDBA(REFORMAT.out)
    MEGAHIT(IN_fastq_raw, megahitKmerSize)
    METAHIPMER2(REFORMAT.out, metahipmer2KmerSize)
    METASPADES(IN_fastq_raw, metaspadesKmerSize)
    MINIA(IN_fastq_raw, miniaKmerSize)
    SKESA(IN_fastq_raw)
    SPADES(IN_fastq_raw, spadesKmerSize)
    UNICYCLER(IN_fastq_raw)
    VELVETOPTIMISER(IN_fastq_raw)

    emit:
    all_assemblies = ABYSS.out.assembly | mix(BCALM2.out.assembly, 
                                            GATBMINIAPIPELINE.out.assembly,
                                            IDBA.out.assembly,
                                            MEGAHIT.out.assembly,
                                            METAHIPMER2.out.assembly,
                                            METASPADES.out.assembly,
                                            MINIA.out.assembly,
                                            SKESA.out.assembly,
                                            SPADES.out.assembly,
                                            UNICYCLER.out.assembly,
                                            VELVETOPTIMISER.out.assembly)
    all_versions = ABYSS.out.version | mix(BCALM2.out.version, 
                                        GATBMINIAPIPELINE.out.version,
                                        IDBA.out.version,
                                        MEGAHIT.out.version,
                                        METAHIPMER2.out.version,
                                        METASPADES.out.version,
                                        MINIA.out.version,
                                        SKESA.out.version,
                                        SPADES.out.version,
                                        UNICYCLER.out.version,
                                        VELVETOPTIMISER.out.version) | collect

}

