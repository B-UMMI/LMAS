# LMAS postprocessing module

This module contains the following processes:

- PROCESS_COMPLETNESS
  - Computes plot with breadth of coverage by number of aligned contigs for each reference
- PLOT_LX 
  -  Computes plot of L metric, from 0 to 100%, each assembly by each reference
- PLOT_NAX
  -  Computes plot of NA metric, from 0 to 100%, each assembly by each reference
- PLOT_NGX
  - Computes plot of NG metric, from 0 to 100%, each assembly by each reference
- PROCESS_SHRIMP_PLOT
  - Computes the plot of the Pls metric for each contig by each reference
- PLOT_CONTIG_DISTRIBUTION
  - Computes plot of contig size distribution
- GAP_ASSESSMENT
  - Computes the coordinates of the reference not covered by the assembly after mapping
- PLOT_GAP_BOXPLOT
  - Computes plot of the gap lenght distribution
- PLOT_GAP_REFERENCE
  - Computes the plot of the location of sequence not covered by the assembly in the reference
- SNP_ASSESSMENT
  - Computes the location of SNPs in the reference sequence in reference to the assembled sequence
- PLOT_SNP_REFERENCE
  - Computes the plot of the location of SNPs by the assembly in the reference
- MISASSEMBLY
  - Computes the misassembly detection algorithm from the alignment paf files. 
- PROCESS_MISASSEMBLY
  - Computes the json with number of misassembled contigs and misassembly events
- PLOT_MISASSEMBLY
  - Computes plot of misassembled contigs 

It emits the following:

- completness_json
  - plot, in json format, with breadth of coverage by number of aligned contigs for each reference
- contig_distribution_json
  - plot, in json format, of contig size distribution
- lx_json
  - plot, in json format, of L metric
- phred_json
  - plot, in json format, of pls metric and contig size by conting
- gap_reference_json
  - plot, in json format, of the location of gaps in the reference
- snp_reference_json
  - plot, in json format, of the location of SNPs in the reference
- gap_boxplot_json
  - plot, in json format, of the gap size distribution
- misassembly_json
  - json with number of misassemled contigs and misassembly events per assembly
- misassembly_report_json
  - json with number of misassemled contigs and misassembly events per assembly
- nax_json
  - plot, in json format, of NA metric
- ngx_json
  - plot, in json format, of NG metric
- misassembly_reference_json
  - Json with number of misassembled contings and misassembly events per reference
- misassembly_plot_json
  - plot, in json format, of number of blocks and contig size of misassembled contigs

