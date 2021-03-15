Metrics
=======

Global metrics are presented both for the **Original Assembly** and in the **Filtered Assembly**.
Assemblies are filtered by a minumum contig size, defined in the `minLength` parameter. 

Global Metrics
--------------

Contigs
^^^^^^^
Total number of contigs in the **Original Assembly** and in the **Filtered Assembly**, by minimum contig size.
**Default:** 1000 basepairs

Basepairs
^^^^^^^^^
Total number of basepairs in the **Original Assembly** and in the **Filtered Assembly**, by minimum contig size.
**Default:** 1000 basepairs

Maximum Contig Size
^^^^^^^^^^^^^^^^^^^
Lenght on the largest contig in the **Original Assembly**.

Nx
^^^^^^^
Length for which the collection of all contigs of that length or longer, for **Original Assembly** and **Filtered Assembly**, 
covers at least *x* % of the total length of the assembled contigs.
**Default:** 50

Mapped Reads
^^^^^^^^^^^^
Percentage of mapped reads to the **Original Assembly**.

Number of 'N' Basepairs
^^^^^^^^^^^^^^^^^^^^^^
Total of 'N's in assembly basepairs, for **Original Assembly** and **Filtered Assembly**.

Missassembled Contigs
^^^^^^^^^^^^^^^^^^^^^^
Number of missassembled contigs in the **Filtered Assembly**.

Per Reference Metrics
---------------------

Contiguity
^^^^^^^^^^^
Longest single alignment between the assembly and the reference, relative to the reference length.

Multiplicity
^^^^^^^^^^^^
Ratio of the length of alignable assembled sequence to covered sequence on the reference.

Validity
^^^^^^^^^
Ratio of the length of alignable assembled sequence to total basepairs in the aligned contings.

Parsimony
^^^^^^^^^
Cost of the assembly (multiplicity over validity).

Identity
^^^^^^^^^
Ratio of identical basepairs in all aligned contigs to the reference.

Lowest Identity
^^^^^^^^^^^^^^^
Identity of the lowest scoring contig to the reference.

Breadth of Coverage
^^^^^^^^^^^^^^^^^^^
Ratio of covered sequence on the reference by aligned contigs.

Aligned Contigs
^^^^^^^^^^^^^^^
Number of aligned contigs to the reference.

Missassembled Contigs
^^^^^^^^^^^^^^^^^^^^^^
Number of aligned contigs with misassemblies.

Lx
^^
Minimal number of contigs that cover *x* % of the sequence of the reference.
**Default:** 90 

NAx
^^^
Length for which the collection of all aligned contigs of that length or longer covers at least *x* % of the total 
length of the aligned assembled contigs.
**Default:** 50

NGx
^^^
Length for which the collection of all aligned contigs of that length or longer covers at least *x* % of the sequence 
of the reference.
**Default:** 50

Aligned Basepairs
^^^^^^^^^^^^^^^^^
Total basepairs aligned to to the reference.

Number of 'N' Basepairs
^^^^^^^^^^^^^^^^^^^^^^
Total of 'N's in basepairs aligned to to the reference. 