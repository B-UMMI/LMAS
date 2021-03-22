Installation
============

Before installing LMAS, a few dependencies must be installed in your system.

Nextflow
--------

`Nextflow <https://www.nextflow.io/>`_ (version 20.01.0 or higher) can be used on any POSIX compatible system (Linux, OS X, etc). 
It requires BASH and Java 8 (or higher) to be installed. 

.. important::

    Instructions on how to install Nextflow are available `here <https://www.nextflow.io/docs/latest/getstarted.html>`_

Container engine
----------------

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

Clone LMAS
-----------

You can clone this repository with git.

.. code-block:: bash

    git clone https://github.com/cimendes/LMAS.git 

All files will be in your local machine.
The main execution file for Nextflow is ``LMAS.nf``. 