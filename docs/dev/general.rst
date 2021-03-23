General Orientation
===================

LAMS code is in two repositories: `source repository <https://github.com/cimendes/LMAS>`_ and 
`report repository <https://github.com/cimendes/lmas_report>`_. 

LMAS repository
----------------

This repository contains the code for LMAS workflow. The workflow is developed in `Nextflow <https://www.nextflow.io/>`_.
The **main file** is `LMAS.nf <https://github.com/cimendes/LMAS/blob/main/LMAS.nf>`_, to be executed by Nextflow.

Configuration
::::::::::::::

The **main configuration** file is `nextflow.config <https://github.com/cimendes/LMAS/blob/main/nextflow.config>`_ and contains
the main configuration parameters for the execution of LMAS. The **parameters, containers and resources** configuration files, 
that can be altered by the user to adapt LMAS execution, are 
`params.config <https://github.com/cimendes/LMAS/blob/main/params.config>`_, 
`containers.config <https://github.com/cimendes/LMAS/blob/main/containers.config>`_,
and `resources.config <https://github.com/cimendes/LMAS/blob/main/resources.config>`_, respectively. 

Information on how to adjust these values is available `here <../user>`_.

Templates
::::::::::

The `templates <https://github.com/cimendes/LMAS/tree/main/templates>`_ folder contains the custom python scripts used
by LMAS to process the data and compute the evaluation metrics, and are called in the 
`LMAS.nf <https://github.com/cimendes/LMAS/blob/main/LMAS.nf>`_ file. 


Resources
:::::::::

The `resources <https://github.com/cimendes/LMAS/tree/main/resources>`_ folder contains the compiled source code 
for the LMAS report. The report code is available in the `report repository <https://github.com/cimendes/lmas_report>`_.

Lib
::::

The `lib <https://github.com/cimendes/LMAS/tree/main/lib>`_ folder contains custum Groovy code used by LMAS for 
the ``--help`` function. 

Docker
::::::

The dockerfile for the main LMAS container, including all necessary python dependencies for the custom code in the 
`templates <https://github.com/cimendes/LMAS/tree/main/templates>`_  is available in the 
`docker folder <https://github.com/cimendes/LMAS/tree/main/docker/LMAS>`_ 


LMAS report
-------------

This `repository <https://github.com/cimendes/lmas_report>`_. contains the source code for the interactive report that 
comes pre-packaged with LMAS.
This project uses ``npx webpack`` to compile a standalone ``main.js`` file that is integrated into LMAS.
The necessary dependencies for the project are provided in the ``environment.yml`` file available in this repo.


