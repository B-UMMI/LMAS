Installation
============

LMAS can be installed through Github (https://github.com/cimendes/LMAS).
It requires a `Nextflow <https://www.nextflow.io/>`_ installation (version ≥ 21.04.1) 
and can be used on any POSIX compatible system (Linux, OS X, etc). All components of LMAS are executed in `Docker containers<https://www.docker.com/>`_, 
being a container engine required. 

Nextflow allows integration with multiple alternatives, such as `Shifter <https://github.com/NERSC/shifter/>`_ or 
`Singularity <https://singularity.hpcng.org/>`_, so a particular one isn’t required. 

To ensure the robustness of custom python code for the quality assessment of assemblies, **continuous integration** of the python templates 
is performed with `pytest <https://docs.pytest.org/en/6.2.x/>`_ and `GitHub Actions <https://github.com/features/actions>`_. 

Below it's a step by step guide on how to install LMAS and all it's dependencies.

Step 1. Nextflow
-----------------

`Nextflow <https://www.nextflow.io/>`_ (version 20.01.0 or higher) can be used on any POSIX compatible system (Linux, OS X, etc). 
It requires BASH and Java 8 (or higher) to be installed. 

.. important::

    Instructions on how to install Nextflow are available `here <https://www.nextflow.io/docs/latest/getstarted.html>`_

Step 2. Container engine
-------------------------

All components of LMAS are executed in docker containers, which means that you’ll need to have a container engine 
installed. The container engines available are the ones supported by Nextflow:

- `Docker`_,
- `Singularity`_
- Shifter (undocumented)

If you already have any one of these installed, you are good to go as the provided docker containers are compatible 
with all engines available. If not, you’ll need to install one.


Singularity
:::::::::::

Singularity is available to download and install `here <http://singularity.lbl.gov/install-linux>`_.
Make sure that you have singularity v2.5.x or higher.
Note that singularity should be installed as root and available on the machine(s) that
will be running the nextflow processes.

.. important::

    Singularity is available as a bioconda package. However, conda installs singularity
    in user space without root privileges, which may prevent singularity images from
    being correctly downloaded. **Therefore it is not recommended that you install
    singularity via bioconda**.

Docker
::::::

Docker can be installed following the instructions on the website:
https://www.docker.com/community-edition#/download.
To run docker as anon-root user, you'll need to following the instructions
on the website: https://docs.docker.com/install/linux/linux-postinstall/#manage-docker-as-a-non-root-user

Step 3. Clone LMAS
-------------------

You can clone this repository with git.

.. code-block:: bash

    git clone https://github.com/cimendes/LMAS.git 

All files will be in your local machine.
The main execution file for Nextflow is ``LMAS.nf``. 