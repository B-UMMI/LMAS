# LMAS assembly module

This module contains the following processes:

- REFORMAT
  - merges the paired-end read files with [BBtools reformat.sh](https://sourceforge.net/projects/bbmap/)
- ABYSS 
  -  Assemble the input fastq data with [ABySS](https://github.com/bcgsc/abyss) assembler. Returns it's assembly and the software version.  
- GATBMINIAPIPELINE
  - Assemble the input fastq data with [GATB-Minia Pipeline](https://github.com/GATB/gatb-minia-pipeline) assembler. Returns it's assembly and the software version.
- IDBA
  - Assemble the input fastq data with [IDBA-UD](https://github.com/loneknightpy/idba) assembler. Returns it's assembly and the software version.
- MEGAHIT
  - Assemble the input fastq data with [MEGAHIT](https://github.com/voutcn/megahit) assembler. Returns it's assembly and the software version.
- METAHIPMER2
  - Assemble the input fastq data with [MetaHipMer2](https://bitbucket.org/berkeleylab/mhm2) assembler. Returns it's assembly and the software version.
- METASPADES
  - Assemble the input fastq data with [metaSPAdes](https://github.com/ablab/spades) assembler. Returns it's assembly and the software version.
- MINIA
  - Assemble the input fastq data with [MINIA](https://github.com/GATB/minia) assembler. Returns it's assembly and the software version.
- SKESA
  - Assemble the input fastq data with [SKESA](https://github.com/ncbi/SKESA) assembler. Returns it's assembly and the software version.
- SPADES
  - Assemble the input fastq data with [SPAdes](https://github.com/ablab/spades) assembler. Returns it's assembly and the software version.
- UNICYCLER
  - Assemble the input fastq data with [Unicycler](https://github.com/rrwick/Unicycler) assembler. Returns it's assembly and the software version.
- VELVETOPTIMISER
  - Assemble the input fastq data with [VelvetOptimiser](https://github.com/tseemann/VelvetOptimiser) assembler. Returns it's assembly and the software version.


It emits the following:

- all_assemblies
  - Mixed channel containing the sample name, the assembler name and the resulting assembly file (fasta) for all assemblers in the workflow
- all_versions
  - Mixed channel containing the file with assembler version for all assemblers in the workflow
