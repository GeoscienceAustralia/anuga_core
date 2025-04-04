
.. _install_anuga:

Install ANUGA
-------------

Once you have a working `python` environment you can install a prebuilt 
version of ANUGA from `conda-forge` (if you have an Anaconda or Miniforge install) 
or via `pip` if you have a standard python install. 

.. note::
    If you want the most recent update of ANUGA or intend to develop ANUGA code you 
    should install ANUGA from source 
    (see `Install ANUGA for Developers <install_anuga_developers.html>`_).
    



Anaconda or Miniforge
~~~~~~~~~~~~~~~~~~~~~

It is always recommended that you create a separate `conda` environment for 
your ANUGA installation. 

So create a python 3.12 conda environment called `anuga_env` 
(or what ever name you like):

.. code-block:: bash

    conda create -c conda-forge -n anuga_env python=3.12 anuga mpi4py
    conda activate anuga_env

Note we have also installed `mpi4py` to allow anuga to run in parallel. 
On some systems you may need to manually install `mpi4py` to match the 
version of `mpi` you are using.


This has setup and activated a `conda` environment `anuga_env` which is using python 3.12. 
(This version of ANUGA has be tested on 3.9 - 3.13)    

We are now ready to use ANUGA. 



Test ANUGA
~~~~~~~~~~  

You can test your ANUGA installation by running the unit tests via:

.. code-block:: bash

    pytest --pyargs anuga


.. note::

    You will nedd to `activate` the `anuga_env` environment each time you want to use ANUGA.
    
    If you are using standard python you use the `source anuga_env/bin/activate` command to 
    activate the environment.

    If you are using `conda` you use the `conda activate anuga_env` command to activate 
    the environment.

    You can add these activate command to your `.bashrc` or `.bash_profile` file to
    automatically activate the environment when you open a terminal.
