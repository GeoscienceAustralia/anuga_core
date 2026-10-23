#!/bin/bash
#
# Install ANUGA with GPU offloading support using the NVIDIA HPC SDK (nvc).
#
# Run this script as follows:
#
#   bash /path/to/anuga_core/tools/install_anuga_nvc.sh
#
# The script locates nvc in the NVIDIA HPC SDK, activates the conda
# environment, and builds the GPU-enabled anuga package.
#
# Environment variables (all optional):
#
#   PY          Python version to use (default: 3.12, matching
#               install_miniforge.sh).  Ignored when a conda environment is
#               already activated -- that one is used.
#   GPU_ARCH    GPU compute capability, e.g. cc120 (RTX 5070 Blackwell),
#               cc86 (RTX 30xx Ampere), cc90 (H100), cc80 (A100), cc70 (V100).
#               Default: DETECTED from the GPU in this machine via nvidia-smi,
#               falling back to a multi-architecture build if that fails.
#               Getting this wrong matters: the build now contains only the
#               architectures asked for, so a mismatched value produces a
#               package that compiles fine and then crashes on every kernel
#               launch.
#   GPU_AWARE_MPI  1 = build with -Dgpu_aware_mpi=true (default 0).
#               Halo buffers are then allocated with omp_target_alloc and staged
#               with omp_target_memcpy, instead of mapped host buffers moved
#               with 'target update'.  NOTE this does NOT pass device pointers
#               to MPI: anuga/shallow_water/gpu/gpu_halo.c stages every halo
#               through host memory in both paths, because some UCX transports
#               (uct_mm intra-node shared memory) SIGSEGV on device pointers.
#               So it needs no CUDA-aware MPI, and gives no GPUDirect.
#   SKIP_TESTS  1 = do not run the GPU test suite at all.
#   FORCE_TESTS 1 = run it even when no GPU is visible (default: skip, so an
#               HPC login node does not spend a long time on tests that
#               cannot pass).
#   NVHPC_ROOT  Override path to NVIDIA HPC SDK root if auto-detection fails.
#               e.g. /opt/nvidia/hpc_sdk/Linux_x86_64/26.3
#
# Example — build for an A100 with Python 3.13:
#
#   PY=3.13 GPU_ARCH=cc80 bash /path/to/anuga_core/tools/install_anuga_nvc.sh
#
# Prerequisites:
#   - NVIDIA HPC SDK installed (see KNOWN_ISSUES.md for apt install recipe)
#   - conda environment anuga_env_${PY} created via install_miniforge.sh

# Keep this default in step with tools/install_miniforge.sh: the documented flow
# is to run that script and then this one, and if the two disagree this one
# looks for an environment the other never created.
PY=${PY:-"3.12"}

# GPU_ARCH is resolved AFTER nvc is located (see below): choosing it needs to
# know which CUDA toolkits this HPC SDK actually ships, not just which GPU is
# present.  "auto" means "work it out".
GPU_ARCH=${GPU_ARCH:-auto}

set -e

# Keep a transcript -- this is the file to ask for when someone reports that a
# GPU build "installed fine but every test crashes".
LOGFILE=${LOGFILE:-"$HOME/anuga_gpu_install_$(date +%Y%m%d_%H%M%S).log"}
exec > >(tee -a "$LOGFILE") 2>&1
echo "# Logging this installation to: $LOGFILE"

trap 'echo ""; echo "#====================================================="; echo "# Installation failed at line $LINENO"; echo "#====================================================="; exit 1' ERR

SCRIPT=$(realpath "$0")
SCRIPTPATH=$(dirname "$SCRIPT")
ANUGA_CORE_PATH=$(realpath "$SCRIPTPATH/..")

echo "#============================================================"
echo "# ANUGA NVC GPU build  (PY=${PY}  GPU_ARCH=${GPU_ARCH})"
echo "#============================================================"
echo " "

# ------------------------------------------------------------------
# Locate nvc
# ------------------------------------------------------------------
if [ -n "$NVHPC_ROOT" ]; then
    NVC="$NVHPC_ROOT/compilers/bin/nvc"
elif command -v nvc >/dev/null 2>&1; then
    # An nvc already on PATH is the user's choice -- on HPC systems it comes
    # from a module and is often a wrapper (e.g. Gadi:
    # /apps/nvidia-hpc-sdk/wrappers/nvc), not the /opt SDK layout below.
    NVC=$(command -v nvc)
else
    # Search /opt/nvidia/hpc_sdk/Linux_x86_64/ for the newest installed version
    NVHPC_BASE="/opt/nvidia/hpc_sdk/Linux_x86_64"
    if [ -d "$NVHPC_BASE" ]; then
        # Pick the highest version directory
        NVHPC_VER=$(ls "$NVHPC_BASE" | grep -E '^[0-9]+\.[0-9]+$' | sort -V | tail -1)
        NVC="$NVHPC_BASE/$NVHPC_VER/compilers/bin/nvc"
    else
        NVC=""
    fi
fi

if [ -z "$NVC" ] || [ ! -x "$NVC" ]; then
    echo "#=====================================================";
    echo "# ERROR: nvc not found."
    echo "#"
    echo "# On an HPC system, load it as a module first, e.g.:"
    echo "#"
    echo "#   module avail nvhpc          # see what is available"
    echo "#   module load nvhpc"
    echo "#"
    echo "# Otherwise install the NVIDIA HPC SDK:"
    echo "#"
    echo "#   curl -fsSL https://developer.download.nvidia.com/hpc-sdk/ubuntu/DEB-GPG-KEY-NVIDIA-HPC-SDK \\"
    echo "#     | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg"
    echo "#   echo 'deb [signed-by=/usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg] \\"
    echo "#     https://developer.download.nvidia.com/hpc-sdk/ubuntu/amd64 /' \\"
    echo "#     | sudo tee /etc/apt/sources.list.d/nvhpc.list"
    echo "#   sudo apt-get update -y && sudo apt-get install -y nvhpc"
    echo "#"
    echo "# Or set NVHPC_ROOT=/path/to/hpc_sdk/Linux_x86_64/<version>"
    echo "#=====================================================";
    exit 1
fi

echo "# nvc found: $NVC"
"$NVC" --version

# The C++ compiler must come from the same SDK as the C one.
#
# meson.build declares 'cpp' as a project language, so meson always probes for
# a C++ compiler -- and with only CC set it picks the system g++. That is
# harmless today because no C++ source is compiled, but meson.build adds
# '-mp=gpu,multicore' as a LINK argument for language 'cpp' as well as 'c'
# (the nvidia_hpc branch), and g++ does not understand that flag. The first C++
# source added to the tree would therefore fail to link, with an error pointing
# at a flag nobody set by hand.
#
# nvc++ lives beside nvc in every SDK layout, including the module wrappers.
if [ -z "$CXX" ]; then
    NVCXX="$(dirname "$NVC")/nvc++"
    if [ -x "$NVCXX" ]; then
        export CXX="$NVCXX"
        echo "# nvc++ found: $CXX"
    else
        echo "# WARNING: nvc++ not found beside nvc; meson will fall back to the"
        echo "#          system C++ compiler. Harmless while the tree has no C++"
        echo "#          sources, but see the note above if that changes."
    fi
else
    echo "# CXX already set, leaving it alone: $CXX"
fi
echo " "

# ------------------------------------------------------------------
# Resolve GPU_ARCH=auto, then PROVE the choice compiles before spending
# several minutes on the real build.
#
# Two things constrain the answer, and neither is guessable:
#   * which GPU is present (nvidia-smi) -- absent on a login node;
#   * which CUDA toolkits this SDK ships.  cc70 (V100) needs CUDA <= 12.x;
#     an SDK carrying only CUDA 13 rejects it outright ("A CUDA toolkit
#     matching the current driver version ... was not installed"), because
#     CUDA 13's libnvvm dropped Volta.
# ------------------------------------------------------------------
arch_compiles() {
    # Tiny OpenMP-target probe; ~2s. Returns 0 if nvc accepts this -gpu= string.
    local arch="$1" tmp rc
    tmp=$(mktemp -d)
    printf '#include <stdio.h>\nint main(void){double a[10];\n#pragma omp target teams distribute parallel for map(tofrom:a[0:10])\nfor(int i=0;i<10;i++)a[i]=i;\nprintf("%%f\\n",a[9]);return 0;}\n' > "$tmp/probe.c"
    "$NVC" -mp=gpu -gpu="$arch" -c "$tmp/probe.c" -o "$tmp/probe.o" >"$tmp/log" 2>&1
    rc=$?
    [ $rc -ne 0 ] && PROBE_ERR=$(head -1 "$tmp/log")
    rm -rf "$tmp"
    return $rc
}

# nvidia-smi is the only way we have to see a GPU, and it is not always cheap.
# Being ON PATH says nothing: an HPC login node typically ships the binary
# without a driver for it to talk to, where it fails with
#
#     NVIDIA-SMI has failed because it couldn't communicate with the NVIDIA
#     driver.
#
# but takes its time about it. The script used to call nvidia-smi three times
# with no timeout, right after announcing the test suite, which is
# indistinguishable from a hung test run. Probe ONCE, bounded, and reuse it.
nvidia_smi_query() {
    command -v nvidia-smi >/dev/null 2>&1 || return 1
    if command -v timeout >/dev/null 2>&1; then
        timeout 10 nvidia-smi --query-gpu="$1" --format=csv,noheader 2>/dev/null
    else
        nvidia-smi --query-gpu="$1" --format=csv,noheader 2>/dev/null
    fi
}

if command -v nvidia-smi >/dev/null 2>&1; then
    echo "# Looking for a GPU (nvidia-smi)..."
fi
# Require BOTH a zero exit status and a plausible device name. Exit status
# alone is not enough (some builds report success with an empty list), and a
# non-empty line alone is not enough either: a failing nvidia-smi may still
# print something -- "[Not Supported]", "[N/A]", or an error -- on stdout,
# which a bare -n test reads as a GPU and sends the script on to run the GPU
# test suite on a machine that has none.
# NB capture the status BEFORE piping: `$?` after a pipeline is the status of
# the last element, so `nvidia_smi_query ... | head -1` would always report
# head's 0 and the check below would never fire.
GPU_PROBE_RAW="$(nvidia_smi_query name)"
GPU_PROBE_RC=$?
GPU_NAME_PROBE="$(printf '%s\n' "$GPU_PROBE_RAW" | head -1)"
case "$GPU_NAME_PROBE" in
    ''|'['*|*'Not Supported'*|*'N/A'*|*'failed'*|*'Error'*|*'error'*)
        GPU_NAME_PROBE="" ;;
esac
[ "$GPU_PROBE_RC" -ne 0 ] && GPU_NAME_PROBE=""
if [ -n "$GPU_NAME_PROBE" ]; then
    HAVE_GPU=1
    echo "# GPU detected: ${GPU_NAME_PROBE}"
else
    HAVE_GPU=0
    if command -v nvidia-smi >/dev/null 2>&1; then
        echo "# No usable GPU here: nvidia-smi is installed but reports no"
        echo "# device (no driver loaded, or none allocated to this node)."
        echo "# This is normal on an HPC login node. Building anyway."
    else
        echo "# No GPU here: nvidia-smi is not installed. Building anyway."
    fi
fi

GPU_ARCH_EXPLICIT=1
[ "$GPU_ARCH" = "auto" ] && GPU_ARCH_EXPLICIT=0

if [ "$GPU_ARCH" = "auto" ]; then
    # Only ask again if the first probe actually saw a device; otherwise this
    # is a second guaranteed-to-fail call, and on a driverless node a second
    # wait for it.
    if [ "$HAVE_GPU" = "1" ]; then
        CAP=$(nvidia_smi_query compute_cap | head -1 | tr -d ' ')
    else
        CAP=""
    fi
    if [[ "$CAP" =~ ^[0-9]+\.[0-9]+$ ]]; then
        GPU_ARCH="cc${CAP//./}"
        echo "# GPU_ARCH=auto -> ${GPU_ARCH} (detected: compute capability ${CAP})"
    else
        # No GPU visible: an HPC login node, or a driver that is not loaded.
        # Build for everything this SDK can target, so the result runs on the
        # compute nodes.
        #
        # Which architectures those are cannot be read off the filesystem: nvc
        # may be a module wrapper (Gadi: /apps/nvidia-hpc-sdk/wrappers/nvc) with
        # no SDK layout beneath it.  cc70 (V100) additionally needs a CUDA <=
        # 12.x, since CUDA 13's libnvvm dropped Volta -- and when an SDK ships
        # both, the CUDA 13 default must be overridden explicitly.  So ASK the
        # compiler: try the widest list first and keep the first that builds.
        echo "# GPU_ARCH=auto -> no GPU visible here; probing what this nvc can build."
        echo "#"
        echo "#   NOTE: this is the SLOW path, and it is the usual one on a login"
        echo "#   node. Probing compiles a test file once per candidate list, and"
        echo "#   the build that follows then generates device code for EVERY"
        echo "#   architecture in the winning list -- up to seven. Expect tens of"
        echo "#   minutes, with long silences from nvc. It is working, not hung."
        echo "#"
        echo "#   If you know the compute nodes you are targeting, naming them is"
        echo "#   much faster, e.g. GPU_ARCH=cc80 (A100) or GPU_ARCH=cc90 (H100)."
        echo "#"
        ALL="cc75,cc80,cc86,cc89,cc90,cc120"
        GPU_ARCH=""
        for CANDIDATE in "cc70,${ALL}" \
                         "cuda12.9,cc70,${ALL}" \
                         "cuda12.8,cc70,${ALL}" \
                         "cuda12.6,cc70,${ALL}" \
                         "${ALL}"; do
            printf '#   trying %-34s ... ' "$CANDIDATE"
            if arch_compiles "$CANDIDATE"; then
                echo "ok"
                GPU_ARCH="$CANDIDATE"
                break
            fi
            echo "no"
        done
        if [ -z "$GPU_ARCH" ]; then
            echo "#====================================================="
            echo "# ERROR: nvc rejected every architecture list tried."
            echo "#        Last error: ${PROBE_ERR}"
            echo "#        Set GPU_ARCH explicitly, e.g. GPU_ARCH=cc80"
            echo "#====================================================="
            exit 1
        fi
        case "$GPU_ARCH" in
            *cc70*) echo "#   V100 (cc70) IS included - this build covers V100 through Blackwell." ;;
            *)      echo "#   NOTE: cc70/V100 could not be built with this SDK (CUDA 13 dropped"
                    echo "#         Volta). This build will NOT run on a V100; load an nvhpc"
                    echo "#         module with CUDA 12.x if you need one." ;;
        esac
    fi
fi

echo "# Checking that nvc accepts -gpu=${GPU_ARCH} ..."
if ! arch_compiles "$GPU_ARCH"; then
    echo "#   rejected: ${PROBE_ERR}"
    # Dropping an architecture the user asked for by name would hand them a
    # build that crashes on the very GPU they were targeting, so only the
    # auto-resolved list is narrowed. An explicit GPU_ARCH is their decision to
    # correct.
    GPU_ARCH_NO70=$(echo "$GPU_ARCH" | sed 's/cuda12\.[0-9]*,//; s/cc70,//; s/,cc70//')
    if [ "$GPU_ARCH_EXPLICIT" = "1" ] && [ "$GPU_ARCH_NO70" != "$GPU_ARCH" ]; then
        echo "#====================================================="
        echo "# ERROR: you asked for ${GPU_ARCH}, but this nvc cannot build cc70."
        echo "#        ${PROBE_ERR}"
        echo "#"
        echo "# cc70 (V100) needs a CUDA <= 12.x toolkit; CUDA 13 dropped Volta."
        echo "# Either load an nvhpc module built against CUDA 12.x, or drop"
        echo "# cc70 yourself and accept that the build will not run on a V100:"
        echo "#"
        echo "#     GPU_ARCH=${GPU_ARCH_NO70} bash ${SCRIPT}"
        echo "#====================================================="
        exit 1
    fi
    if [ "$GPU_ARCH_NO70" != "$GPU_ARCH" ] && arch_compiles "$GPU_ARCH_NO70"; then
        echo "#   retrying without cc70 (V100): ${GPU_ARCH_NO70}"
        echo "#   NOTE: the result will NOT run on a V100."
        GPU_ARCH="$GPU_ARCH_NO70"
    else
        echo "#====================================================="
        echo "# ERROR: nvc cannot build for -gpu=${GPU_ARCH}"
        echo "#        ${PROBE_ERR}"
        echo "#"
        echo "# Pick architectures this SDK supports, e.g."
        echo "#     GPU_ARCH=cc80,cc90 bash ${SCRIPT}"
        echo "#====================================================="
        exit 1
    fi
fi
echo "#   ok - building for ${GPU_ARCH}"
echo " "

# ------------------------------------------------------------------
# Select the conda environment
#   Prefer an already-activated environment (this is what the user is working
#   in — likely with anuga already installed). Otherwise fall back to a
#   miniforge install in $HOME with env anuga_env_$PY.
# ------------------------------------------------------------------
if [ -n "$CONDA_PREFIX" ] && [ -n "$CONDA_DEFAULT_ENV" ] \
        && [ "$CONDA_DEFAULT_ENV" != "base" ]; then
    ENV_NAME="$CONDA_DEFAULT_ENV"
    CONDA_RUN=""   # build/test run in the current, already-activated shell
    echo "# Using the active conda environment: ${ENV_NAME}  (${CONDA_PREFIX})"
else
    CONDA_BIN="$HOME/miniforge3/bin"
    if [ ! -f "$CONDA_BIN/conda" ]; then
        echo "#=====================================================";
        echo "# ERROR: no conda environment is active and miniforge3 was not"
        echo "#        found at $HOME/miniforge3."
        echo "#"
        echo "# Either activate your anuga environment first"
        echo "#   (e.g. conda activate anuga_env_${PY}), or run"
        echo "#   install_miniforge.sh to create one."
        echo "#=====================================================";
        exit 1
    fi
    ENV_NAME="anuga_env_${PY}"
    if ! "$CONDA_BIN/conda" env list | grep -q "${ENV_NAME}"; then
        FOUND=$("$CONDA_BIN/conda" env list | awk '/^anuga_env_/ {print $1}' | tr '\n' ' ')
        echo "#=====================================================";
        echo "# ERROR: conda environment '${ENV_NAME}' not found."
        if [ -n "$FOUND" ]; then
            echo "#"
            echo "# These anuga environments do exist: ${FOUND}"
            echo "# Activate one, or name it explicitly, e.g."
            echo "#     PY=<version> bash ${SCRIPT}"
        else
            echo "# Run install_miniforge.sh first (PY=${PY}), or activate your env."
        fi
        echo "#=====================================================";
        exit 1
    fi
    CONDA_RUN="$CONDA_BIN/conda run -n ${ENV_NAME}"
    echo "# conda environment: ${ENV_NAME}  (via $HOME/miniforge3)"
fi
echo " "

# ------------------------------------------------------------------
# Preflight: build backend + build requirements
#
# The build below uses `pip install --no-build-isolation`, so pip does NOT
# create a temporary build environment — the meson-python backend (module
# `mesonpy`) and the rest of pyproject's build-system.requires must already be
# installed in the target environment.  If they are missing, pip dies with an
# opaque "BackendUnavailable: Cannot import 'mesonpy'" traceback that never
# mentions meson-python.  Check up front and say exactly what to install.
#
# This also reports the environment's *actual* Python version: when an env is
# already activated we use it and ignore $PY, so the banner's PY can differ.
# ------------------------------------------------------------------
echo "# Preflight: checking build backend and build requirements"

PREFLIGHT_PY='
import importlib.util, shutil, sys
missing = []
for mod, pkg in (("mesonpy", "meson-python"), ("Cython", "cython"),
                 ("pybind11", "pybind11"), ("numpy", "numpy")):
    if importlib.util.find_spec(mod) is None:
        missing.append(pkg)
for exe in ("meson", "ninja"):
    if shutil.which(exe) is None:
        missing.append(exe)
print("PREFLIGHT|%d.%d|%s" % (sys.version_info[0], sys.version_info[1],
                              " ".join(missing)))
'
PREFLIGHT_OUT=$($CONDA_RUN python -c "$PREFLIGHT_PY" 2>/dev/null || true)
PREFLIGHT_LINE=$(printf '%s\n' "$PREFLIGHT_OUT" | grep '^PREFLIGHT|' | tail -1 || true)

if [ -z "$PREFLIGHT_LINE" ]; then
    echo "#====================================================="
    echo "# ERROR: could not run python in environment '${ENV_NAME}'."
    echo "#        Is the environment usable?  Try:  ${CONDA_RUN} python -V"
    echo "#====================================================="
    exit 1
fi

ENV_PY_VER=$(printf '%s' "$PREFLIGHT_LINE" | cut -d'|' -f2)
MISSING=$(printf '%s' "$PREFLIGHT_LINE" | cut -d'|' -f3)

echo "#   environment '${ENV_NAME}' is Python ${ENV_PY_VER}"

if [ -n "$MISSING" ]; then
    echo "#====================================================="
    echo "# ERROR: environment '${ENV_NAME}' is missing build requirements:"
    echo "#"
    echo "#     ${MISSING}"
    echo "#"
    echo "# This build uses 'pip install --no-build-isolation', so the"
    echo "# meson-python backend and its build requirements must already be"
    echo "# installed in the environment.  Without them pip fails with an"
    echo "# opaque \"BackendUnavailable: Cannot import 'mesonpy'\"."
    echo "#"
    echo "# Fix - install them into this environment:"
    echo "#"
    echo "#   conda install -c conda-forge ${MISSING}"
    echo "#"
    echo "# Or recreate the environment with everything already in it:"
    echo "#"
    echo "#   conda env create -n anuga_env_${ENV_PY_VER} \\"
    echo "#       -f environments/environment_${ENV_PY_VER}.yml"
    echo "#   conda activate anuga_env_${ENV_PY_VER}"
    echo "#====================================================="
    exit 1
fi

echo "#   ok - meson-python, meson, ninja, cython, pybind11, numpy all present"
echo " "

# ------------------------------------------------------------------
# Build ANUGA with GPU offloading
# ------------------------------------------------------------------
echo "#============================================================"
echo "# Building ANUGA with GPU offloading"
echo "#   CC=$NVC"
echo "#   CXX=${CXX:-<system default; no C++ sources today>}"
echo "#   gpu_offload=true  gpu_arch=${GPU_ARCH}"
echo "#   gpu_aware_mpi=${GPU_AWARE_MPI:-0}  (halo BUFFER allocation only --"
echo "#     ANUGA never hands device pointers to MPI, on either setting, so"
echo "#     this says nothing about whether your MPI is CUDA-aware.  A 0 here"
echo "#     is not a failed detection: it is the default, and it is unrelated"
echo "#     to what 'ompi_info | grep cuda' reports.  See issue #223.)"
echo "#============================================================"
echo " "

cd "${ANUGA_CORE_PATH}"

# meson-python reuses the build/cp<ver> directory and only honours CC on the
# FIRST configure of a dir — a subsequent build just runs `meson setup
# --reconfigure`, which keeps the originally detected compiler. A leftover dir
# configured with gcc (e.g. a prior CPU build) therefore stays on gcc and the
# gpu_offload=true guard in meson.build rejects it ("not supported with gcc").
# Remove any stale build dir so CC=nvc takes effect on a clean configure.
echo "# Removing any stale meson build directory (build/cp*) for a clean nvc configure"
rm -rf "${ANUGA_CORE_PATH}"/build/cp*
echo " "

# EDITABLE=0 installs a copy into site-packages instead, which survives having
# the source tree or its build directory removed.  The default stays editable
# (this script is mostly used on a development checkout), but note the
# dependency: an editable install imports from ${ANUGA_CORE_PATH} and loads its
# compiled extensions from build/cp<ver>.  Delete that directory later and
# `import anuga` fails with a FileNotFoundError naming a missing build path,
# which does not obviously mean "reinstall".
if [ "${EDITABLE:-1}" = "1" ]; then
    PIP_TARGET_ARGS="-e ."
    echo "# Installing EDITABLE (in place). Keep ${ANUGA_CORE_PATH}/build/cp* -"
    echo "# removing it breaks 'import anuga'. Use EDITABLE=0 for a standalone copy."
else
    PIP_TARGET_ARGS="."
    echo "# Installing a COPY into site-packages (EDITABLE=0)."
fi
echo " "

GPU_AWARE_MPI_ARG=""
if [ "${GPU_AWARE_MPI:-0}" = "1" ]; then
    GPU_AWARE_MPI_ARG="-Csetup-args=-Dgpu_aware_mpi=true"
    echo "# GPU_AWARE_MPI=1: building with -Dgpu_aware_mpi=true"
    echo "#   (device-allocated halo buffers; MPI itself still receives host"
    echo "#    pointers - see gpu_halo.c. No CUDA-aware MPI is required.)"
    if ! $CONDA_RUN python -c "import mpi4py" >/dev/null 2>&1; then
        echo "#   WARNING: mpi4py is not installed in this environment, so the"
        echo "#            extension will be built against the single-process"
        echo "#            stubs and this flag will have no effect. Install"
        echo "#            mpi4py (and pymetis) first if you want parallel runs."
    fi
    echo " "
fi

CXX_ARG=""
[ -n "$CXX" ] && CXX_ARG="CXX='$CXX'"

$CONDA_RUN bash -c \
    "CC='$NVC' ${CXX_ARG} pip install --no-build-isolation -v ${PIP_TARGET_ARGS} \
     -Csetup-args=-Dgpu_offload=true \
     -Csetup-args=-Dgpu_arch=${GPU_ARCH} \
     ${GPU_AWARE_MPI_ARG}"

echo " "
# ---------------------------------------------------------------------------
# Report what was built, and stop if it cannot run on this machine.
#
# A build that targets only architectures NEWER than the GPU present compiles
# cleanly and then fails at every kernel launch -- which surfaces as a wall of
# CRASH from the tests below, with nothing explaining why.  (The reverse is
# fine: nvc embeds PTX, which the driver JIT-compiles forward, so a build for an
# older architecture runs on a newer GPU.)
# ---------------------------------------------------------------------------

# ------------------------------------------------------------------
# Post-build check: does this toolchain track header dependencies?
#
# On NVHPC, meson's PGI compiler mixin does not implement
# get_dependency_gen_args(), so the compile line never gets -MD/-MF even though
# meson still writes "deps = gcc" and a depfile into the ninja rule. Ninja then
# records ZERO dependencies for every object and NO header change ever triggers
# a rebuild.
#
# That is a silent-wrong-binary bug, not a slow-build bug. struct gpu_domain
# embeds struct domain, so when only some translation units are rebuilt they
# disagree about member offsets: observed as fixed_flux_timestep aliasing
# evolve_max_timestep, which made mode 2 print "Using a fixed timestep!" and
# hang, with no compile error. Editing a header does it; so does merely
# switching to a branch whose headers differ, which is easier to hit because
# nobody edited anything.
#
# This script always does a clean configure, so THIS build is correct either
# way. The check exists to say whether a plain incremental
# `pip install -e .` is safe afterwards.
# ------------------------------------------------------------------
BUILD_DIR=$(ls -d "${ANUGA_CORE_PATH}"/build/cp* 2>/dev/null | head -1)
if [ -n "$BUILD_DIR" ] && command -v ninja >/dev/null 2>&1; then
    DEP_TOTAL=$(ninja -C "$BUILD_DIR" -t deps 2>/dev/null \
                | grep -cE '^[^ ]+\.o: #deps [1-9]' || true)
    if [ "${DEP_TOTAL:-0}" -eq 0 ]; then
        echo "#=================================================================="
        echo "# WARNING: this build has NO header dependency tracking."
        echo "#"
        echo "#   ninja recorded zero header dependencies for every object, so a"
        echo "#   changed header will NOT trigger a rebuild. An incremental"
        echo "#   'pip install --no-build-isolation -e .' can then produce a"
        echo "#   SILENTLY WRONG binary whose translation units disagree about"
        echo "#   struct layout - no compile error, no crash, wrong answers or"
        echo "#   a hang."
        echo "#"
        echo "#   Cause: meson's PGI/NVHPC compiler mixin does not implement"
        echo "#   get_dependency_gen_args(), so -MD/-MF is never passed. nvc"
        echo "#   supports the flags; meson simply does not ask for them."
        echo "#"
        echo "#   UNTIL THAT IS FIXED UPSTREAM: after ANY header edit or branch"
        echo "#   switch, re-run this script. Never a plain incremental"
        echo "#   'pip install -e .'."
        echo "#=================================================================="
        echo " "
    else
        echo "# Header dependency tracking: OK (${DEP_TOTAL} objects with recorded deps)."
        echo "#   Incremental 'pip install -e .' is safe after a header change."
        echo " "
    fi
fi

echo "#============================================================"
echo "# Build report"
echo "#============================================================"
if ! $CONDA_RUN python "${ANUGA_CORE_PATH}/tools/anuga_build_report.py" --check; then
    echo ""
    echo "#====================================================="
    echo "# Stopping before the tests: the build above cannot run"
    echo "# on this machine, so every GPU test would report CRASH."
    echo "#====================================================="
    exit 1
fi
echo " "

# Do not run the GPU tests where there is no GPU.  The usual case is an HPC
# login node: every test would fail or skip after a long wait, and running a
# heavy suite there is antisocial (many sites forbid it outright).  The build is
# still complete and usable -- the tests just have to happen where the GPUs are.
if [ "${SKIP_TESTS:-0}" = "1" ]; then
    echo "SKIP_TESTS=1 - skipping the GPU test suite."
elif [ "$HAVE_GPU" = "0" ] && [ "${FORCE_TESTS:-0}" != "1" ]; then
    echo "#=================================================================="
    echo "# Skipping the GPU test suite: no GPU is visible here."
    echo "#"
    echo "# The build itself is complete. These tests need a GPU, and this"
    echo "# looks like a login node - running them here would take a long"
    echo "# time, fail anyway, and load a shared machine."
    echo "#"
    echo "# Run them where the GPUs are, in a job or interactive session:"
    echo "#"
    echo "#   python ${ANUGA_CORE_PATH}/scripts/anuga_run_isolated_tests.py"
    echo "#"
    echo "# and check the build suits that node's GPU with:"
    echo "#"
    echo "#   python ${ANUGA_CORE_PATH}/tools/anuga_build_report.py --check"
    echo "#"
    echo "# FORCE_TESTS=1 runs them here regardless."
    echo "#=================================================================="
else
    echo "#============================================================"
    echo "# Running GPU test suite (isolated runner)"
    echo "#   One fresh process per test.  A plain 'pytest' on this file"
    echo "#   auto-skips on a GPU build: the NVHPC OpenMP-target runtime"
    echo "#   aborts once many mode-2 GPU domains are created in a single"
    echo "#   process, so the tests must each run in their own process."
    echo "#   scripts/anuga_run_isolated_tests.py defaults to test_DE_gpu_omp.py"
    echo "#   and opts in via ANUGA_GPU_TESTS_ISOLATED=1.  Run the script"
    echo "#   directly (not the installed console command, which an editable"
    echo "#   'pip install -e .' does not place on PATH)."
    echo "#============================================================"
    echo " "
    $CONDA_RUN \
        python "${ANUGA_CORE_PATH}/scripts/anuga_run_isolated_tests.py"
fi

echo " "
echo "#=================================================================="
echo "# Congratulations! ANUGA GPU build succeeded."
echo "#"
echo "# To use GPU mode, activate the environment and select the 'unified'"
echo "# compute mode (mode 2).  Either per-domain in Python:"
echo "#"
echo "#   conda activate ${ENV_NAME}"
echo "#   python -c \\"
echo "#     \"import anuga; d = anuga.rectangular_cross_domain(100,100); \\"
echo "#      d.set_boundary({b: anuga.Reflective_boundary(d) for b in d.get_boundary_tags()}); \\"
echo "#      d.set_compute_mode('unified')\""
echo "#"
echo "# or process-wide (applies to every domain) via the environment:"
echo "#"
echo "#   export ANUGA_DEFAULT_COMPUTE_MODE=unified   # 'legacy' (CPU) | 'unified' (CPU/GPU)"
echo "#"
echo "# On this GPU build, 'unified' offloads to the GPU -- confirm with"
echo "# 'python -c \"import anuga; print(anuga.gpu_offload_enabled())\"' (True)."
echo "# 'legacy' runs on the CPU.  (d.set_multiprocessor_mode(2) is the old"
echo "# alias for 'unified' and still works.)"
echo "#"
echo "# Under MPI, run one rank per GPU: 'unified' with more ranks than GPUs"
echo "# oversubscribes the device and deadlocks."
echo "#=================================================================="
