Help with porting to Python3:


#--------------
# Web resources
#--------------
https://docs.python.org/3/howto/pyporting.html
http://python3porting.com/

# On Unicode strings, str, bytes, bytearray
https://learning-python.com/strings30.html

#-------------
# Git commands
#-------------


# Main repository is: git@github.com:stoiver/anuga_core.git

# Synchronise porting branch from upstream
git remote add upstream git@github.com:stoiver/anuga_core.git
git pull upstream python3_anuga_py3 

#---------
# Burndown
#---------

Date PR Result
-------------------

20200527     Ran 1035 tests in 39.265s FAILED (errors=463, failures=42)
20200528  37 Ran 1035 tests in 39.265s FAILED (errors=272, failures=42)
20200528  38 Ran 1034 tests in 37.107s FAILED (errors=254, failures=43)
20200529  39 Ran 1146 tests in 79.962s FAILED (errors=234, failures=40) Hashing, checksums and exec(s, globals())
20200529     Ran 1146 tests in 76.251s FAILED (errors=233, failures=40) char_to_string (netcdf)
20200529     Ran 1146 tests in 73.962s FAILED (errors=243, failures=28) mesh boundaries text vs binary
