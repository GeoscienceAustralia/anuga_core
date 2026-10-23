#!/bin/bash
#
# Run this script as follows:
#
# bash /path/to/anuga_core/tools/install_miniforge.sh
#
# This script will install miniforge3 and create a conda environment called
# anuga_env_${PY} where PY is the python version you want to use.
# The script will then install the anuga package from the anuga_core directory
# and run the unittests.
#
# By default a python version of 3.12 will be installed. 
#
# If you want to install a different version of python, set the PY environment
# variable before running the script. For example, to install python 3.9 run the
# script as follows:
#
# PY=3.9 bash /path/to/anuga_core/tools/install_miniforge.sh
#
# Then the script will install python 3.9 and create the anuga_env_3.9 environment.
#
# This builds the CPU (gcc) package, which is the right one for CPU and MPI
# runs. A GPU build is a SEPARATE build with different trade-offs: it offloads
# to the GPU, but its CPU fallback is single-threaded and much slower than this
# one (see docs/source/appendices/install_gpu.rst). The script therefore does
# NOT switch compilers just because a GPU is present -- it tells you the GPU
# path exists and leaves the choice to you.
#
# To do both in one go, add GPU=1:
#
#   GPU=1 bash /path/to/anuga_core/tools/install_miniforge.sh
#
# which runs the CPU install, then chains into tools/install_anuga_nvc.sh.


PY=${PY:-"3.12"}

set -e

# Keep a transcript.  An install that goes wrong is usually reported second-hand
# ("it failed"), and this is the file to ask for.
LOGFILE=${LOGFILE:-"$HOME/anuga_install_$(date +%Y%m%d_%H%M%S).log"}
exec > >(tee -a "$LOGFILE") 2>&1
echo "# Logging this installation to: $LOGFILE"


trap 'echo ""; echo "#====================================================="; echo "# Installation failed at line $LINENO"; echo "#====================================================="; exit 1' ERR


SCRIPT=$(realpath "$0")
SCRIPTPATH=$(dirname "$SCRIPT")
ANUGA_CORE_PATH=$(realpath "$SCRIPTPATH/..")


#test PY>3.8 and <3.14
if [[ "$PY" =~ ^3\.(1[0-4]|[9])$ ]]; then
     echo "Requested python version is $PY"
     echo " "
else
    echo "Python version must be greater/equal than 3.9 and less than 3.15"
    exit 1
fi


echo "#==========================="
echo "# Install miniforge3"
echo "#==========================="
cd $HOME

if [ -d "$HOME/miniforge3" ]; then
     echo "miniforge3 seems to already exist."
else
     echo "miniforge3 does not exist."
     echo "Installing from script Miniforge3.sh..."
     if [ -f "$HOME/Miniforge3.sh" ]; then
          echo "Running Miniforge3.sh first..."
     else
          echo "Miniforge3.sh does not exist. Downloading..."
          wget -O "$HOME/Miniforge3.sh" "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
     fi
     bash Miniforge3.sh -b
fi

echo " "
echo "#==============================================="
echo "# create conda environment anuga_env_${PY}"
echo "#==============================================="
echo "..."
ENV_YML="${SCRIPTPATH}/../environments/environment_${PY}.yml"
if [ ! -f "$ENV_YML" ]; then
     echo "ERROR: no environment file for Python ${PY}: $ENV_YML"
     echo "Available: $(ls ${SCRIPTPATH}/../environments/environment_3.*.yml \
                        | grep -v intel | sed 's#.*environment_##;s#\.yml##' | tr '\n' ' ')"
     exit 1
fi
# Re-running the script must not fail on an environment that already exists:
# `conda env create` errors out in that case, so update it in place instead.
if ./miniforge3/bin/conda env list | grep -qE "^anuga_env_${PY}\s"; then
     echo "Environment anuga_env_${PY} already exists - updating it in place."
     ./miniforge3/bin/conda env update --name anuga_env_${PY} --file "$ENV_YML" --prune
else
     ./miniforge3/bin/conda env create --file "$ENV_YML"
fi

echo " "
echo "#======================================"
echo "# activate environment anuga_env_${PY}"
echo "#======================================"
echo "..."
source ./miniforge3/bin/activate anuga_env_${PY}


echo "#================================================================"
echo "# Installing anuga from the ${ANUGA_CORE_PATH} directory"
echo "#================================================================"
echo "..."

cd ${SCRIPTPATH}
cd ..
# Non-editable (a copy into site-packages) on purpose: it keeps working no
# matter what happens to this source tree or its build/ directory.  Developers
# who want to edit the sources in place can pass EDITABLE=1, but then the
# build/cp<ver> directory must be kept -- deleting it breaks `import anuga`.
if [ "${EDITABLE:-0}" = "1" ]; then
     echo "# EDITABLE=1: installing in place (keep the build/cp* directory!)"
     pip install --no-build-isolation -e .
else
     pip install --no-build-isolation .
fi
echo " "

echo "#==========================="
echo "# Build report"
echo "#==========================="
python "${ANUGA_CORE_PATH}/tools/anuga_build_report.py" || true
echo " "

echo "#==========================="
echo "# Run unittests"
echo "#==========================="
echo "#   Set SKIP_TESTS=1 to skip, or FAST_TESTS=1 for the quick subset."
echo " "

cd sandpit
if [ "${SKIP_TESTS:-0}" = "1" ]; then
     echo "SKIP_TESTS=1 - skipping the test suite."
elif [ "${FAST_TESTS:-0}" = "1" ]; then
     pytest --pyargs anuga --run-fast
else
     pytest --pyargs anuga
fi

echo " "
echo "#=================================================================="
echo "# Congratulations, Looks like you have successfully installed anuga"
echo "#=================================================================="

echo "#=================================================================="
echo "# To use anuga you must activate the python environment anuga_env_${PY}"
echo "# that has just been created. Run the command"
echo "# "
echo "# source ~/miniforge3/bin/activate anuga_env_${PY}"
echo "# "
echo "#=================================================================="

# ---------------------------------------------------------------------------
# GPU: report what is possible on this machine, and act only if asked.
#
# Deliberately NOT automatic. A gpu_offload=true (nvc) build is not a superset
# of this one: NVHPC's host fallback does not scale with threads, so a GPU
# build used on the CPU is several times slower than this gcc build (measured;
# see claude/KNOWN_ISSUES.md). Choosing the compiler from the hardware would
# silently downgrade anyone whose work is CPU or MPI bound, and would also
# break the test step above, which must be run per-process on a GPU build.
# ---------------------------------------------------------------------------
GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)
GPU_CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -d ' ')
NVC_PATH=$(command -v nvc 2>/dev/null)
if [ -z "$NVC_PATH" ]; then
     NVC_PATH=$(ls -d /opt/nvidia/hpc_sdk/Linux_x86_64/*/compilers/bin/nvc 2>/dev/null | sort -V | tail -1)
fi

if [ -n "$GPU_NAME" ]; then
     echo " "
     echo "#=================================================================="
     echo "# A CUDA GPU was detected:"
     echo "#     ${GPU_NAME} (compute capability ${GPU_CAP:-unknown} -> cc${GPU_CAP//./})"
     echo "#"
     echo "# What you have just installed is the CPU (gcc) build. Keep it for"
     echo "# CPU and MPI runs: it is several times faster on the host than a"
     echo "# GPU build is, because NVHPC's host fallback is single-threaded."
     if [ -n "$NVC_PATH" ]; then
          echo "#"
          echo "# nvc is installed (${NVC_PATH}), so you can also build for the GPU:"
          echo "#"
          echo "#     bash ${SCRIPTPATH}/install_anuga_nvc.sh"
          echo "#"
          echo "# It detects cc${GPU_CAP//./} automatically. Note it REPLACES this"
          echo "# build in this environment - the two cannot coexist in one env."
     else
          echo "#"
          echo "# For GPU offload you would also need the NVIDIA HPC SDK (nvc),"
          echo "# which is not installed. See docs/source/appendices/install_gpu.rst."
     fi
     echo "#=================================================================="
fi

if [ "${GPU:-0}" = "1" ]; then
     echo " "
     echo "#=================================================================="
     echo "# GPU=1 - continuing with the GPU (nvc) build"
     echo "#=================================================================="
     if [ -z "$NVC_PATH" ]; then
          echo "#====================================================="
          echo "# ERROR: GPU=1 was requested but nvc was not found."
          echo "#        Install the NVIDIA HPC SDK first, or drop GPU=1."
          echo "#        See docs/source/appendices/install_gpu.rst"
          echo "#====================================================="
          exit 1
     fi
     if [ -z "$GPU_NAME" ]; then
          echo "# NOTE: no GPU detected by nvidia-smi. Building anyway (GPU=1);"
          echo "#       install_anuga_nvc.sh will build for every supported"
          echo "#       architecture when it cannot see a GPU."
     fi
     # PY targets the environment just created (though the GPU script will
     # prefer the one this script activated, which is the same one).  LOGFILE is
     # shared so a GPU=1 run leaves a single transcript rather than two.
     PY="${PY}" LOGFILE="${LOGFILE}" bash "${SCRIPTPATH}/install_anuga_nvc.sh"
fi

echo "#=================================================================="
echo "# NOTE: If you run the command"
echo "# "
echo "# conda init"
echo "# "
echo "# (which will change your .bashrc file) "
echo "# then in new terminals you will be able to use "
echo "# the conda command"
echo "# "
echo "# conda activate anuga_env_${PY}"
echo "# "
echo "# to activate the anuga_env_${PY} environment"
echo "#=================================================================="
