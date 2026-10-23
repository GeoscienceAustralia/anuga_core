#!/bin/bash
#
# Install ANUGA and its build dependencies into a Python venv on Ubuntu.
#
# Detects the Ubuntu LTS version and applies version-appropriate numpy/gdal
# pins, then installs the apt build deps, creates a venv, pip-installs the
# Python stack, builds ANUGA, and runs the test suite.
#
#   bash tools/install_ubuntu.sh
#
# For the most robust install (esp. the fiddly numpy/gdal combination), prefer
# the conda path instead: tools/install_miniforge.sh + `pip install -e .`.
# See tools/README.md.
#
# License: 3-clause BSD

set -e

echo "#==========================="
echo "# Determine Ubuntu version"
echo "#==========================="
VERSION_ID=$(grep -oP 'VERSION_ID="\K[\d.]+' /etc/os-release)
echo "Ubuntu $VERSION_ID"
echo " "

# Version-specific pip pins. numpy/gdal are the parts that vary between
# releases; everything else is shared below.
case "$VERSION_ID" in
    22.04)
        NUMPY_PIN="numpy==1.26"
        GDAL_PIN="gdal==3.4.1"
        ;;
    24.04)
        NUMPY_PIN="numpy==2.2"
        GDAL_PIN="gdal"
        ;;
    26.04)
        # Provisional: 26.04 ships Python 3.14, so numpy/gdal are left to pip to
        # resolve. Untested — if the venv build fails, use install_miniforge.sh.
        NUMPY_PIN="numpy"
        GDAL_PIN="gdal"
        ;;
    *)
        echo "Ubuntu $VERSION_ID is not supported by this script."
        echo "Supported: 22.04, 24.04, 26.04. For anything else use the conda"
        echo "path: tools/install_miniforge.sh + 'pip install -e .'."
        exit 1
        ;;
esac

echo "#==========================="
echo "# Install packages via apt"
echo "#==========================="
sudo apt install -y -q build-essential python-dev-is-python3 gfortran netcdf-bin \
    libnetcdf-dev libhdf5-serial-dev gdal-bin libgdal-dev libopenmpi-dev openmpi-bin \
    python3-venv

echo "#======================================="
echo "# Create a virtual environment and then"
echo "# install python packages via pip"
echo "#======================================="
cd "$(dirname "${BASH_SOURCE[0]}")"/..
ANUGA_CORE_PATH=$(pwd)
echo "ANUGA_CORE_PATH: $ANUGA_CORE_PATH"

python3 -m venv anuga_env
source anuga_env/bin/activate

# Ensure meson picks the pip-installed numpy, not system numpy. Derive the
# venv's python minor version so this works across Ubuntu releases (rather than
# hard-coding python3.10 / python3.12 / ...).
PYVER=$(python -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')
export PKG_CONFIG_PATH="${ANUGA_CORE_PATH}/anuga_env/lib/python${PYVER}/site-packages/numpy/_core/lib/pkgconfig:${PKG_CONFIG_PATH}"

pip install wheel "$NUMPY_PIN" scipy "$GDAL_PIN" matplotlib pytest cython netcdf4 \
    dill future gitpython pyproj pymetis pybind11 meshpy Pmw ipython \
    utm affine mpi4py xarray meson meson-python ninja

echo "#============================================"
echo "# Install anuga from the anuga_core directory"
echo "#============================================"
pip install .

echo "#==========================="
echo "# Run unittests"
echo "#==========================="
# On numpy 2.x + apt gdal some tests may fail due to a numpy/gdal ABI
# incompatibility; that does not abort the install. For a fully working
# environment use tools/install_miniforge.sh.
pytest -q --pyargs anuga || echo "(some tests failed — see the numpy/gdal note above)"

echo "#================================================"
echo "# Done. Activate the environment with:"
echo "#   source anuga_env/bin/activate"
echo "#================================================"
