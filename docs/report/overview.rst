Overview
========

LMAS creates an interactive HTML report, stored in the "report/" folder in the directory where the workflow was executed. 
To open the report simply click oh the "index.html" file and the report will open on your default browser.

On the top right corner there's a direct link to LMAS `source repository <https://github.com/cimendes/LMAS>`_ and 
`documentation <https://lmas.readthedocs.io/en/latest>`_. 

.. image:: ../resources/report_lmas.gif
    :align: center
    :scale: 70 %


Report Overview
----------------

The top portion of the report contains information on the input samples and overall performance of the assemblers in LMAS.
This is devided into three tabs:
* Overview
* Performance
* About us

.. image:: ../resources/overview.gif
    :align: center
    :scale: 70 %


Overview
::::::::

This is devided into two sections: **Input Data** and **About**. 

Input data presents the name of the file passed to '--reference' parameter, the sample names of the short-read paired-end 
raw sequencing files passed to '--fastq', and the number of reads in each file. 

The About loads information about the samples used, in markdown, passed on to LMAS with the '--md' parameter. This is an 
optional parameter so, if missing, no information will be presented. 


Performance
:::::::::::

This tab has a table with information on version, container used and performance metrics for each assembler in LMAS.

* **Assembler:** Assembler name
* **Version:** Version of the assembler used, obtained through the '--version' or '-v' flag. If no version is provided, the cell is left blank.
* **Container:** Container used in LMAS for the assemblers.
* **Avg Time:** Average exectution time to assemble the input samples. 
* **CPU/Hour:** Average amount of time the CPU was used by the Assembler 
* **Max Memory (GB):** Maximum amount of memory, in Gigabytes, used by the Assembler.
* **Average Read (GB):** Average data read, in Gigabytes, by the Assembler.
* **Average Write (GB):** Average data written, in Gigabytes, by the Assembler.

About Us
::::::::

Contains information of LMAS development team. 

