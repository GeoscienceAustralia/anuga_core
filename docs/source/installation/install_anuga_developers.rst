
Install ANUGA for Developers
----------------------------

If you want to use the very latest version of ANUGA (or develop ANUGA code) then you need
to download the `anuga_core` repository from `github` and then `pip` install 
ANUGA from the source. These steps will require that the following package `git` is installed.


First download the `anuga_core` repository from `github` and then run our `install_miniforge.sh`
script.

Download ANUGA from `github`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We need to download (clone) the ANUGA source code from `github`

.. code-block:: bash

    git clone https://github.com/anuga-community/anuga_core.git

This creates a directory `anuga_core`.

.. note::

    If you want to also contribute to the code base, you must have a GitHub 
    account and setup authentication from your developer workstation to GitHub 
    as per these instructions:  https://docs.github.com/en/authentication/managing-commit-signature-verification. The command to clone ANUGA as a developer is then 

    .. code-block:: bash

        git clone git@github.com:anuga-community/anuga_core.git

Install ANUGA using Script
~~~~~~~~~~~~~~~~~~~~~~~~~~~

We have a script in the `anuga_core/tools` directory that will install `Miniforge` 
and ANUGA and its dependencies.

Simply run the following command from the `anuga_core` directory:

.. code-block:: bash

    bash tools/install_miniforge.sh

This will create a `conda` python 3.12 environment `anuga_env_3.12` and install ANUGA 
and its dependencies.

.. note::

    If you want to install ANUGA for a different version of python, you can set the PY 
    environment variable when running the `install_miniforge.sh` as follows:
    
    
    .. code-block:: bash

      export PY=3.11; bash tools/install_miniforge.sh
    
    This will install ANUGA for python 3.11. 

.. note::

    The script `install_miniforge.sh` essentially does the following:

    .. code-block:: bash

        cd anuga_core
        conda env create -n anuga_env_3.12 -f environments/environment.3.12.yml
        conda activate anuga_env_3.12

    and finally installs ANUGA in editable mode via: 

    .. code-block:: bash

        pip install --no-build-isolation -e .


.. note::

    You may need to install a compiler to complete the `pip install`. 
    You can use the system compilers or use `conda` to install compilers as such:

    .. code-block:: bash

        conda install compilers

    or for win32:

    .. code-block:: bash

        conda install m2w64-gcc libpython 

    or for macOS:

    For macOS we suggest installing `homebrew` which will 
    provide the gcc compilers. Once you have installed
    the compilers via `homebrew` you need to set the environment variables to
    point to the `gcc` and `g++` compilers. 
    
    Check their location via:
    
    .. code-block:: bash
        
        which gcc
        which g++

    and then set the environment variables as such:

    .. code-block:: bash

        export CC=/opt/homebrew/bin/gcc
        export CXX=/opt/homebrew/bin/g++

    or whatever the path is to your homebrew compilers.

    Once you have installed the compilers you can run the `pip install` command
    to install ANUGA.

    .. code-block:: bash

        pip install --no-build-isolation -e .


Testing the installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Once the installation is complete you can activate the `anuga_env_3.12` environment
and run the unit tests to check that everything is working. 

Test the installation.

.. code-block:: bash

   pytest --pyargs anuga


Updating
~~~~~~~~

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

  pytest -q --pyargs anuga 

 


