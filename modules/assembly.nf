process ABYSS {

    tag { sample_id }
    label 'process_assembly'
    publishDir "results/$sample_id/assembly/abyss/"

    when:
    params.abyss

    input:
    tuple sample_id, file(fastq) from IN_ABYSS
    val KmerSize from Channel.value(params.abyssKmerSize)
    val BloomSize from Channel.value(params.abyssBloomSize)

    output:
    tuple sample_id, val('ABySS'), file('*_ABySS.fasta') into OUT_ABYSS
    file '.*version' into ABYSS_VERSION

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
    tuple sample_id, file(fastq) from IN_BCALM2
    val KmerSize from Channel.value(params.bcalmKmerSize)

    output:
    tuple sample_id, val('BCALM2'), file('*_BCALM2.fasta') into OUT_BCALM2
    file '.*version' into BCALM2_VERSION

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
    tuple sample_id, file(fastq_pair) from IN_GATB_MINIA_PIPELINE
    val kmer_list from Channel.value(params.gatbKmerSize)
    val do_error_correction from GATB_error_correction
    val besst_iter from Channel.value(params.gatb_besst_iter)

    output:
    tuple sample_id, val('GATBMiniaPipeline'), file('*_GATBMiniaPipeline.fasta') into OUT_GATB
    file '.*version' into GATB_VERSION

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

process reformat_IDBA {

    tag { sample_id }
    label 'process_assembly'

    when:
    params.idba

    input:
    tuple sample_id, file(fastq_pair) from IN_IDBA

    output:
    tuple sample_id, file('*.fasta') into REFORMAT_IDBA

    script:
    "reformat.sh in=${fastq_pair[0]} in2=${fastq_pair[1]} out=${sample_id}_reads.fasta"
}

process IDBA {

    tag { sample_id }
    label 'process_assembly'
    publishDir "results/$sample_id/assembly/IDBA-UD/"

    when:
    params.idba

    input:
    tuple sample_id, file(fasta_reads_single) from  REFORMAT_IDBA

    output:
    tuple sample_id, val('IDBA-UD'), file('*_IDBA-UD.fasta') into OUT_IDBA
    file '.*version' into IDBA_VERSION

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
    publishDir "results/$sample_id/assembly/MEGAHIT/", pattern: '*_megahit*.fasta'

    when:
    params.megahit

    input:
    tuple sample_id, file(fastq_pair) from IN_MEGAHIT
    val kmers from Channel.value(params.megahitKmerSize)

    output:
    tuple sample_id, val('MEGAHIT'), file('*_MEGAHIT.fasta') into OUT_MEGAHIT
    file '.*version' into MEGAHIT_VERSION

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

process reformat_METAHIPMER2 {

    tag { sample_id }
    label 'process_assembly'

    when:
    params.metahipmer2

    input:
    tuple sample_id, file(fastq_pair) from IN_METAHIPMER2

    output:
    tuple sample_id, file('*.fastq') into REFORMAT_METAHIPMER2

    script:
    "reformat.sh in=${fastq_pair[0]} in2=${fastq_pair[1]} out=${sample_id}_reads.fastq"
}

process METAHIPMER2 {

    tag { sample_id }
    label 'process_assembly'
    publishDir "results/$sample_id/assembly/MetaHipMer2/"

    when:
    params.metahipmer2

    input:
    tuple sample_id, file(fasta_reads_single) from  REFORMAT_METAHIPMER2
    val kmer from Channel.value(params.metahipmer2KmerSize)

    output:
    tuple sample_id, val('MetaHipMer2'), file('*_MetaHipMer2.fasta') into OUT_METAHIPMER2
    file '.*version' into METAHIPMER2_VERSION

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
    tuple sample_id, file(fastq_pair) from IN_METASPADES
    val kmers from Channel.value(params.metaspadesKmerSize)

    output:
    tuple sample_id, val('metaSPAdes'), file('*_metaspades.fasta') into OUT_METASPADES
    file '.*version' into METASPADES_VERSION

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
    tuple sample_id, file(fastq) from IN_MINIA
    val kmer from Channel.value(params.miniaKmerSize)

    output:
    tuple sample_id, val('MINIA'), file('*_minia.fasta') into OUT_MINIA
    file '.*version' into MINIA_VERSION

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
    tuple sample_id, file(fastq_pair) from IN_SKESA

    output:
    tuple sample_id, val('SKESA'), file('*_skesa.fasta') into OUT_SKESA
    file '.*version' into SKESA_VERSION

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
    tuple sample_id, file(fastq_pair) from IN_SPADES
    val kmers from Channel.value(params.spadesKmerSize)

    output:
    tuple sample_id, val('SPAdes'), file('*_spades.fasta') into OUT_SPADES
    file '.*version' into SPADES_VERSION

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
    tuple sample_id, file(fastq_pair) from IN_UNICYCLER

    output:
    tuple sample_id, val('Unicycler'), file('*_unicycler.fasta') into OUT_UNICYCLER
    file '.*version' into UNICYCLER_VERSION

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
    publishDir "results/$sample_id/assembly/VelvetOtimiser"

    when:
    params.velvetoptimiser

    input:
    tuple sample_id, file(fastq_pair) from IN_VELVETOPTIMISER

    output:
    tuple sample_id, val('VelvetOptimiser'), file('*.fasta') into OUT_VELVETOPTIMISER
    file '.*version' into VELVETOPTIMISER_VERSION

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