
Install ANUGA for Developers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you want to use the very latest version of ANUGA (or to develop ANUGA code) then you need
to download the `anuga_core` repository from `github` and then `pip` install
ANUGA from the source. These steps will require that the following package `git` is installed.

ANUGA supports **Python 3.10 – 3.14**.  Environment files for each supported
version are provided under ``environments/`` (e.g. ``environment_3.12.yml``).


The process involves downloading the `anuga_core` repository from `github` and then running the `install_miniforge.sh` 
(or `install_miniforge_windows.bat` on windows) script from the `anuga_core/tools` directory. 
This will install `Miniforge` (if not already installed) and then create a `conda` environment with the required dependencies 
and finally `pip` install ANUGA in editable mode via the `-e` option of the `pip install` command.

Here are the details.

Download ANUGA from `github`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We need to download (clone) the ANUGA source code from `github`

.. code-block:: bash

    git clone https://github.com/anuga-community/anuga_core.git

This creates a directory `anuga_core`.

.. note::

    If you want to also contribute to the code base, you must have a GitHub 
    account and setup authentication from your developer workstation to GitHub 
    as per these instructions:  https://docs.github.com/en/authentication/managing-commit-signature-verification. 
    The command to clone ANUGA as a developer is then 

    .. code-block:: bash

        git clone git@github.com:anuga-community/anuga_core.git

Choose a branch
^^^^^^^^^^^^^^^

The clone leaves you on ``main``, which is the released line -- currently
ANUGA 4.0.0. Work that has not been released yet lives on ``develop``:

.. code-block:: bash

    cd anuga_core
    git checkout develop

Take ``develop`` if you want the passive tracer and sediment transport
modules, :ref:`tracers` and :ref:`sediment`, which are not in 4.0.0 and are
the parts we would most like tested before the next release. Stay on ``main``
if you want the released code.

.. note::

    The documentation is published for both branches, and they describe
    different code:

    * https://anuga.readthedocs.io/en/latest/ follows ``main``
    * https://anuga.readthedocs.io/en/develop/ follows ``develop``

    Read the one that matches the branch you built, or the API you find in
    the documentation may not be the API you have installed.

Install ANUGA using Script
^^^^^^^^^^^^^^^^^^^^^^^^^^

We have a scripts in the `anuga_core/tools` directory that will install `Miniforge` 
and ANUGA and its dependencies.

Simply run the following command from the `anuga_core` directory:

.. code-block:: bash

    bash tools/install_miniforge.sh

or on windows run the script:

.. code-block:: bash

    call tools\\install_miniforge_windows.bat

This will create a `conda` python 3.12 environment `anuga_env_3.12` and install ANUGA 
and its dependencies.

.. note::

    If you want to install ANUGA for a different version of python, you can set the PY 
    environment variable when running the `install_miniforge.sh` as follows:
    
    
    .. code-block:: bash

      export PY=3.11; bash tools/install_miniforge.sh

    or for windows:

    .. code-block:: bash

      set PY=3.11 && call tools\\install_miniforge_windows.bat
    
    This will install ANUGA for python 3.11. 

.. note::

    The install scripts essentially does the following:

    1. Downloads and Installs `Miniforge` if not already installed.
    2. Creates a `conda` environment with the required dependencies for ANUGA (including required compilers on linux, windows and macos).
    3. Installs ANUGA in editable mode via `pip install --no-build-isolation -e .`

    For instance for python 3.12 the script  essentially does the following (without error checking):

    .. code-block:: bash

      wget -O "$HOME/Miniforge3.sh" "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
      bash "$HOME/Miniforge3.sh" -b -p "$HOME/miniforge3"
      cd anuga_core
      conda env create -n anuga_env_3.12 -f environments/environment_3.12.yml
      conda activate anuga_env_3.12
      conda install compilers
      pip install --no-build-isolation -e .


.. note::

    A compiler is needed to complete the `pip install`. 
    You can use the system compilers or use `conda` to install compilers as such:

    For linux:

    .. code-block:: bash

        conda install compilers

    or for win32:

    .. code-block:: bash

        conda install libpython gcc_win-64 gxx_win-64

    or for macOS:

    .. code-block:: bash

        conda install cxx-compiler llvm-openmp

    Once you have installed the compilers you can run the `pip install` command
    to install ANUGA.

    .. code-block:: bash

        pip install --no-build-isolation -e .

    The `--no-build-isolation` option is needed to ensure that the dependencies (in particular the compilers)
    installed in the `conda` environment are used during the build process.

MPI parallel support
^^^^^^^^^^^^^^^^^^^^

To run ANUGA simulations in parallel across multiple processes (using
``mpiexec`` / ``mpirun``), two additional packages are required:

* **mpi4py** — Python bindings for MPI.  Provides the communication layer
  used by ``anuga.distribute``, ``anuga.distribute_basic_mesh``, and related
  functions.
* **pymetis** — Python wrapper for the METIS graph-partitioning library.
  Used by ANUGA to decompose the mesh into balanced subdomains before
  distributing to MPI ranks.

Both packages are already included in the conda environment files under
``environments/``.  For example, ``environment_3.12.yml`` lists them at
lines 16 and 24 respectively, so the standard ``conda env create`` step
installs them automatically.

If you need to add them to an existing environment:

.. code-block:: bash

    conda activate anuga_env_3.12
    conda install -c conda-forge mpi4py pymetis

Verify MPI is working:

.. code-block:: bash

    python -c "from mpi4py import MPI; print('MPI size:', MPI.COMM_WORLD.Get_size())"
    mpiexec -np 4 python -c "from mpi4py import MPI; print(MPI.COMM_WORLD.Get_rank())"

See :doc:`../parallel/index` for how to write and run parallel ANUGA scripts.

.. note::

    On HPC clusters, load the system MPI module *before* activating the
    conda environment so that ``mpi4py`` links against the optimised
    system MPI rather than the bundled conda-forge one::

        module load openmpi/4.1.5   # use your cluster's module name
        conda activate anuga_env_3.12


Testing the installation
^^^^^^^^^^^^^^^^^^^^^^^^

Once the installation is complete you can activate the ``anuga_env_3.12`` environment
and run the unit tests to check that everything is working.

.. code-block:: bash

    cd sandpit
    conda activate anuga_env_3.12
    pytest --pyargs anuga

This runs the full test suite (~1 600 tests, approximately 3 minutes).

Fast test run
"""""""""""""

For a quick feedback loop during development, pass ``--run-fast`` to skip the
slow tests (MPI-based parallel tests and a handful of individually marked
long-running tests):

.. code-block:: bash

    pytest --pyargs anuga --run-fast

This runs approximately 1 500 tests in around 40 seconds.

Tests are considered *slow* if they:

* live under ``anuga/parallel/tests/`` (MPI tests that spawn subprocesses), or
* are decorated with ``@pytest.mark.slow`` (individually identified long-running
  tests in other modules).

To run *only* the slow tests:

.. code-block:: bash

    pytest --pyargs anuga -m slow

To mark a new test as slow, add the decorator before the test method::

    import pytest

    @pytest.mark.slow
    def test_my_expensive_computation(self):
        ...

ANUGA also comes with a validation test suite which verifies the correctness of
real life hydraulic scenarios. You can run them as follows:

.. code-block:: bash

    cd validation_tests
    python run_auto_validation_tests.py

Using the installation
^^^^^^^^^^^^^^^^^^^^^^

You can now use ANUGA by activating the `anuga_env_3.12` environment and then running your python scripts
that use ANUGA.

.. code-block:: bash

    conda activate anuga_env_3.12
    python my_anuga_script.py

If you have a machine with multiple cores you might want to run your anuga scripts using multiple threads. For instance 4 threads.
You can set the environment variable `OMP_NUM_THREADS=4`, as such:

.. code-block:: bash

    conda activate anuga_env_3.12
    export OMP_NUM_THREADS=4
    python my_anuga_script.py

.. note::

   Run your scripts from a directory other than the repository root. Python
   searches the working directory before the installed packages, so a script
   started from inside the ``anuga_core`` clone imports the source tree in
   preference to the installed ANUGA. If that tree has not been built (or was
   built for a different environment) the import fails with a message saying
   ``anuga/_version.py`` is missing; the fix is to ``cd`` somewhere else, or
   to build the tree in place with ``pip install --no-build-isolation -e .``.


Updating
^^^^^^^^^

From time to time you might like to update your version of anuga to the latest version on 
github. You can do this by going to the `anuga_core` directory and `pulling` the latest
version and then reinstalling via the following commands:
 
.. code-block:: bash

  conda activate anuga_env_3.12
  cd anuga_core
  git pull
  pip install --no-build-isolation -editable .

And finally check the new installation by running the unit tests via:

.. code-block:: bash
    
    cd sandpit
    pytest -q --pyargs anuga

 


