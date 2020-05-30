=============================
Help with porting to Python3:
=============================

--------------
Web resources
--------------
https://docs.python.org/3/howto/pyporting.html
http://python3porting.com/

# On Unicode strings, str, bytes, bytearray
https://learning-python.com/strings30.html

-------------
 Git commands
-------------

# Main repository is: git@github.com:stoiver/anuga_core.git

# Clone ANUGA
git clone git@github.com:uniomni/anuga_core.git
git checkout python3_anuga_py3    # Py3 branch
git checkout anuga_py3     # Older Py3 branch - works under Python2.7?


# Synchronise porting branch from upstream
git remote add upstream git@github.com:stoiver/anuga_core.git
git pull upstream python3_anuga_py3 

----------------------------------
# Create ANUGA Python environments
conda create -n anuga27 -c conda-forge python=2.7   git pip nose numpy scipy netcdf4 matplotlib gdal dill cython future openmp mpi4py
conda create -n anuga38 -c conda-forge python=3.8.2 git pip nose numpy scipy netcdf4 matplotlib gdal dill cython future openmp mpi4py
----------------------------------



---------
 Burndown
---------


Env: anuga38
------------
Date      PR Result                                                     Notes ::
20200527     Ran 1035 tests in 39.265s FAILED (errors=463, failures=42) 
20200528  37 Ran 1035 tests in 39.265s FAILED (errors=272, failures=42) 
20200528  38 Ran 1034 tests in 37.107s FAILED (errors=254, failures=43)
20200529  39 Ran 1146 tests in 79.962s FAILED (errors=234, failures=40) Hashing, checksums and exec(s, globals())
20200529     Ran 1146 tests in 76.251s FAILED (errors=233, failures=40) char_to_string (netcdf)
20200529  40 Ran 1146 tests in 73.962s FAILED (errors=243, failures=28) mesh boundaries text vs binary
20200529     Ran 1147 tests in 80.033s FAILED (errors=206, failures=12) sorting and removal of duplicates. All 55 mest tests ok.
20200529     Ran 1147 tests in 79.454s FAILED (errors=206, failures=10) All 19 loadASCII tests ok
20200529     Ran 1147 tests in 80.269s FAILED (errors=191, failures=10) Syntax_tree fixed in system tools. All 25 tests ok.
20200531     Ran 1147 tests in 66.287s FAILED (errors=197, failures=10) Replaced assertCountEqual with assertEqual which works in both environments

Env: anuga27
------------
Date      PR Result                                                     Notes ::
20200530     Ran 1290 tests in 96.139s FAILED (errors=13, failures=9)   After work on python3.8.2
20200530     Ran 1290 tests in 103.100s FAILED (errors=5, failures=9)   Fixed some porting bugs in caching
20200531     Ran 1290 tests in 104.523s FAILED (errors=1, failures=9)   Replaced assertCountEqual with assertEqual which works in both environments
