""" Build configuration for ANUGA. """

# Flags to pass to GCC
GCCFLAGS = '-O3 -Wall'

# Flags to pass to MSVC Compiler
MSVCFLAGS = '/Wall'

# Optimisation level used when building the .pyo files
OPTIMISATION_LEVEL = 1

# Where to install the compiled files. Defaults vary based on OS.
# Uncommenting this and leaving it as the empty string will cause the
# build to fail.
#PREFIX = ''

# Install .py files as well as the binaries? Yes/No option
INSTALL_PYTHON_SOURCE = 'Yes'
