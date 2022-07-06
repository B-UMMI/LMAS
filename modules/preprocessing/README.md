# LMAS preprocessing module

This module contains the following processes:

- PROCESS_REFERENCE 
  -  input reference sequence (multifasta allowed) to it's tripled version, so the reference replicon is concatenated 3 times to allow start-to-end overlaps in the downstream mapping processes
- PROCESS_READS
  - Collects information on the number of reads in the input fastq files

It emits the following:

- triple_reference
  - Fasta file with the input reference replicons concatenated 3 times
- reads_info
  - json file with the read number per input files (paired-end fastq)
