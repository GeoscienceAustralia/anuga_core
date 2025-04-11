.. _install_python:

Install Python
==============

ANUGA requires python 3.X (X>8) and a number of other python packages as defined in the 
`requirements.txt <https://github.com/anuga-community/anuga_core/blob/main/requirements.txt>`_ 
file. These packages are installed via `pip` or `conda` as described below.

So the first requirement is to have a working python environment. A clean way 
to install the dependencies necessary for ANUGA is to use the Anaconda, 
or the Miniforge Python distributions by Continuum Analytics.

We recommend either using
`Anaconda` or `Miniforge` python environments. `Miniforge` is a minimal version of `Anaconda`
and is recommended for installing `ANUGA` as it has a smaller footprint.

.. _Install Anaconda:

Install Anaconda
----------------

To install `Anaconda` follow the instructions at
`Anaconda <https://www.anaconda.com/products/individual>`_.

Once you have installed `Anaconda` open a terminal which will open a `base` conda environment.
Goto the section :ref:`Install ANUGA` to create a new conda environment for ANUGA.


.. _Install Miniforge:

Install MiniForge
-----------------

 

Using a `conda` installation has the advantage of allowing you to easily create multiple 
python environments and is particularly 
useful if you want to keep multiple versions of ANUGA

The most stable install is via the `conda-forge` channel
which is easily available using `Miniforge`. In particular the installation of 
the `gdal` and `mpi4py` modules are more stable using `Miniforge`. 
We recommend  using `Miniforge`. 

These conda environments do not require administrative rights 
to your computer and do not interfere with the Python installed in your system. 

So first we should install the latest version of `Miniforge` from  https://github.com/conda-forge/miniforge or
use, for instance, `wget` to download the latest version via:

.. code-block:: bash

    wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"


.. note::
    
   If on ubuntu you don't have `wget` you can install it via: 

   .. code-block:: bash

        sudo apt-get update -q
        sudo apt-get install wget
   
   On MacOS you can install `wget` via `brew` and on Windows you can download the installer from the
   `Miniforge` website.


Run the installation script:

.. code-block:: bash

    bash Miniforge3.sh

and then activate `miniforge` by running 

.. code-block:: bash

    source miniforge3/bin/activate

.. note::

    During the `Miniforge` installation you will be asked to accept the licence 
    (essentially apache 2.0) and whether to run `conda init` to change your `.bashrc` 
    file to allow activation of the 
    base conda environment when opening a new terminal.
    
    If you choose not to run `conda init` you will need to run the 
    following command every time to activate `miniforge`

    .. code-block:: bash

        source miniforge3/bin/activate 


Once `Miniforge` is installed and activated we can now create an environment to run ANUGA. 

