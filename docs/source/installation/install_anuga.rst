
.. _install_anuga:

Install ANUGA
=============

Once you have a working `python` environment you can install a prebuilt 
version of ANUGA from `conda-forge` (if you have an Anaconda or Miniforge install) 
or via `pip` if you have a standard python install. 

.. note::
    If you want the most recent update of ANUGA or intend to develop ANUGA code you
    should install ANUGA from source
    (see :doc:`install_anuga_developers`).
    

Anaconda or Miniforge
---------------------

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


This has set up and activated a ``conda`` environment ``anuga_env`` using Python 3.12.
ANUGA supports Python 3.10 – 3.14.

We are now ready to use ANUGA.


pip
---

If you are using a standard (non-``conda``) Python environment, ANUGA can be
installed from PyPI:

.. code-block:: bash

    pip install anuga

.. note::

    The ``conda-forge`` route above is recommended. ANUGA depends on compiled
    libraries (GDAL, NetCDF, MPI) whose ``pip`` wheels can be harder to get
    working — especially on Windows. If a ``pip install`` cannot build these,
    use the ``conda`` path instead.


Test ANUGA
----------

A quick check that ANUGA imports and reports its version:

.. code-block:: bash

    python -c "import anuga; print(anuga.__version__)"

To run the fast test suite (skips the slow/parallel tests, ~40 s):

.. code-block:: bash

    pytest --pyargs anuga --run-fast

or the full suite (~1600 tests, a few minutes):

.. code-block:: bash

    pytest --pyargs anuga


.. note::

    You will need to `activate` the `anuga_env` environment each time you want to use ANUGA.
    
    If you are using standard python you use the `source anuga_env/bin/activate` command to 
    activate the environment.

    If you are using `conda` you use the `conda activate anuga_env` command to activate 
    the environment.

    You can add the activate command to your `.bashrc` or `.bash_profile` file to
    automatically activate the environment when you open a terminal.
