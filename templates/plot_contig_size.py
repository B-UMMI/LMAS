#!/usr/bin/env python3
"""
Purpose
-------
Script to obtain information of percentage of mapping contigs/basepairs from the filtered assemblies
and produce a descriptive boxplot for mapped and unmapped contigs per assembler.

Expected input
--------------
The following variables are expected whether using NextFlow or the
:py:func:`main` executor.
- ``coverage_files``: list of files with completeness information to be plotted
        e.g.: ``'[SampleA_AssemblerA.csv, SampleB_AssemblerB.csv]'``

Authorship
----------
Rafael Mamede, rmamede@medicina.ulisboa.pt
https://github.com/rfm-targa
Inês Mendes, cimendes@medicina.ulisboa.pt
https://github.com/cimendes
"""

