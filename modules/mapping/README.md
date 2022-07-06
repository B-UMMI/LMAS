# LMAS mapping module

This module contains the following processes:

- FILTER_ASSEMBLY
  - filters out contigs smaller than `--minLength` from an assembly with [BBtools reformat.sh](https://sourceforge.net/projects/bbmap/)
- READ_MAPPING 
  -  Maps the reads to the original and filtered assembly with [minimap2](https://github.com/lh3/minimap2).Returns the percentage of mapped reads for each assembly.
- ASSEMBLY_MAPPING
  -  Maps the filtered assembled contigs to the tripled reference sequences with [minimap2](https://github.com/lh3/minimap2).Returns the paf file in addition the the information recieved as input (sample name, assembly name, filtered assembly).
- ASSEMBLY_STATS_GLOBAL
  - Computes the global statistics for an assembly (number of contigs, number of basepairs, largest contig size, number of uncalled bases, Nx metric and percentage of mapped reads, for original and filtered assembly)
- PROCESS_ASSEMBLY_STATS_GLOBAL
  - Compiles the global statistics into a json file of all assemblies for all samples.
- ASSEMBLY_STATS_MAPPING
  - Computes the reference-dependent statistics for an assembly (lsa, multiplicity, validity, parsimony, identity, lowest identity, breadth of coverage, Lx, number of aligned contigs, NAx, NGx, number of aligned basepairs and number of uncalled bases). Computes the pls for each contig. computes the NA, NG and L metric for an x of 0% to 100%. 
- PROCESS_ASSEMBLY_STATS_MAPPING
  - Compiles the per reference statistics into a json file of all assemblies for all samples.

It emits the following:

- paf
  - file with the mapping information of the filtered assembly to the reference replicon
- stats_global
  - json file with the global statistics
- stats_mapping
  - json file with the per-reference statistics
- boc_csv
  - csv file with number of mapped contigs and breadth of coverage for each assembly by each reference for all samples
- lx_csv
  - csv file with number L metric, from 0 to 100%, each assembly by each reference
- nax_csv
  - csv file with number NA metric, from 0 to 100%, each assembly by each reference 
- ngx_csv
  - csv file with number NG metric, from 0 to 100%, each assembly by each reference 
- phred_csv
  - csv file with the Pls metric for all contigs for each assembly by each reference
- df_csv
  - dataframe with assembly info for each reference by contig

