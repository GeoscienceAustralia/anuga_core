# Meson build file for ANUGA core.
# Version is derived from git tags via `git describe`.
# Tag releases as `3.3.0`, `3.3.1`, etc. — no `v` prefix, matching the
# existing tag convention. Between releases the version will look like
# `3.3.0-4-gabcdef` (4 commits past the tag, short SHA).
#
# Note: project() must be the first statement in meson.build.
# The version is computed inline; git sha/date are fetched after.
project('anuga', 'c', 'cpp', 'cython',
  version: run_command('python', ['_git_version.py'], check: false).stdout().strip(),
  default_options: ['cpp_std=c++17'])

# Fetch additional git metadata for the generated _version.py
git = find_program('git', required: false, native: true)
if git.found()
  _git_sha  = run_command(git, ['rev-parse', 'HEAD'], check: false).stdout().strip()
  _git_date = run_command(git, ['log', '-1', '--format=%ci', 'HEAD'], check: false).stdout().strip()
else
  _git_sha  = 'unknown'
  _git_date = 'unknown'
endif

add_project_arguments('-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION', language: 'cpp')
add_project_arguments('-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION', language: 'c')

py3 = import('python').find_installation(pure: false)

# `meson dist` (which is what meson-python's sdist build runs) has no .git in
# the staged tree, so bake the version in as a file _git_version.py reads
# before trying git. Without this an sdist install reports 0.0.0+unknown.
# Dist scripts run from the build directory, hence files() for the path.
meson.add_dist_script(py3, files('_git_version.py'), '--write-dist')

dep_py = py3.dependency()
py3_version = py3.language_version()

#========================================================
# Setup numpy dependencies
#========================================================
incdir_numpy = run_command(py3,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check: true
).stdout().strip()

if py3_version.version_compare('>3.8')
  # Recommended for Python > 3.8
  dep_numpy = dependency('numpy')
else
  # Fallback for Python 3.8
  inc_numpy = include_directories(incdir_numpy)
  dep_numpy = declare_dependency(include_directories: inc_numpy)
endif



# FIXME SR: Do we use this now?
#message('NumPy include directory: ' + inc_numpy)

dependencies = [dep_py, dep_numpy]

#========================================================
# Deal with OpenMP
#========================================================

cc = meson.get_compiler('c')
openmp = dependency('openmp', required: false)
message('Compiler ID is: ' + cc.get_id())

#========================================================
# Cray Programming Environment Detection
#========================================================
# On Cray systems, cc/CC are wrappers that hide the actual compiler.
# PE_ENV tells us which backend compiler is being used.
# This allows proper flag selection when using Cray compiler wrappers.

pe_env = run_command('sh', ['-c', 'echo $PE_ENV'], check: false).stdout().strip()

if pe_env != ''
  message('Cray Programming Environment detected: ' + pe_env)
  # Map PE_ENV to meson compiler IDs
  if pe_env == 'NVIDIA'
    effective_compiler_id = 'nvidia_hpc'
  elif pe_env == 'GNU'
    effective_compiler_id = 'gcc'
  elif pe_env == 'CRAY'
    effective_compiler_id = 'cray'
  elif pe_env == 'INTEL'
    effective_compiler_id = 'intel'
  elif pe_env == 'AOCC' or pe_env == 'AMD'
    effective_compiler_id = 'clang'
  else
    effective_compiler_id = cc.get_id()
    warning('Unknown PE_ENV value: ' + pe_env + ', falling back to compiler detection')
  endif
  is_cray_wrapper = true
else
  effective_compiler_id = cc.get_id()
  is_cray_wrapper = false
endif

message('Effective compiler ID: ' + effective_compiler_id)

#========================================================
# Quieten nvc (NVHPC) Cython-boilerplate warning noise
#========================================================
# nvc's EDG frontend surfaces warnings that gcc hides with -Wno-unused-function:
# ~1200+ "declared but never referenced" for the __pyx_*/__Pyx_* helpers Cython
# emits into every generated .c, plus the no_sanitize attribute Cython adds that
# nvc doesn't recognise. These are pure boilerplate and bury any real diagnostic.
# get_supported_arguments() silently drops any name this nvc version rejects.
if effective_compiler_id == 'nvidia_hpc'
  add_project_arguments(
    cc.get_supported_arguments([
      '--diag_suppress=declared_but_not_referenced',
      '--diag_suppress=set_but_not_used',
      '--diag_suppress=unrecognized_attribute',
      '--diag_suppress=white_space_inside_splice',
    ]),
    language: 'c')
endif

#========================================================
# Intel MKL Detection
#========================================================
# When using Intel compilers, MKL provides optimised BLAS/LAPACK for numpy/scipy.
# The use_mkl option can be set explicitly; it also auto-enables for Intel compilers.

is_intel_compiler = effective_compiler_id in ['intel', 'intel-cl', 'intel-llvm', 'icx']
use_mkl_opt = get_option('use_mkl')
use_mkl = use_mkl_opt or is_intel_compiler

if use_mkl
  mkl_dep = dependency('mkl-dynamic-lp64-seq', required: false)
  if not mkl_dep.found()
    mkl_dep = dependency('mkl', required: false)
  endif
  if mkl_dep.found()
    dependencies += [mkl_dep]
    message('Intel MKL found and enabled')
  else
    if is_intel_compiler
      warning('Intel compiler detected but MKL not found. ' +
              'Install MKL via: conda install mkl mkl-devel libblas=*=*mkl. ' +
              'Using environment_<version>_intel.yml is recommended.')
    else
      warning('use_mkl=true but MKL not found — install mkl and mkl-devel.')
    endif
  endif
  if is_intel_compiler
    # -qmkl is supported by both classic Intel (icc/icpc) and Intel LLVM (icx/icpx).
    # This flag selects the MKL BLAS/LAPACK implementation at the C-extension level
    # and works alongside the OpenMP flags (-qopenmp / -fiopenmp) set below.
    add_project_arguments('-qmkl=sequential', language: ['c', 'cpp'])
    add_project_link_arguments('-qmkl=sequential', language: ['c', 'cpp'])
    message('MKL compiler flag -qmkl=sequential added for Intel compiler')
  endif
endif


if openmp.found()
  # On Windows, the mingw compiler does not support OpenMP ATOMIC operations
  # so using gcc_win-64 gxx_win-64 provided by conda-forge
  # On linux, OpenMP is supported by gcc and g++
  # On macOS, OpenMP is not supported by the Apple clang compiler so using
  # clang and llvm-openmp provided by conda-forge
  #
  # Note: These flags are for the original mode 1 code (sw_domain_openmp_ext)
  # The new mode 2 code (sw_domain_gpu_ext) has its own flags in shallow_water/meson.build
  #
  # Uses effective_compiler_id to handle Cray wrappers correctly
  #
  # Do NOT use -march=native here: it bakes the build host's ISA into the
  # compiled extensions, so a wheel built on e.g. a Sapphire Rapids CI runner
  # gets AVX-512 and dies with an illegal instruction at `import anuga` on any
  # CPU without it (this broke the 3.3.8 PyPI wheels, e.g. on Colab). Leave the
  # toolchain's portable default in place (conda-forge sets -march=nocona).

  if effective_compiler_id in ['intel', 'intel-cl']
    openmp_c_args = ['-O3', '-qopenmp', '-g']
    add_project_link_arguments('-qopenmp', language: ['c', 'cpp'])
  elif effective_compiler_id in ['intel-llvm', 'icx']
    openmp_c_args = ['-O3', '-fiopenmp', '-g']
    add_project_link_arguments('-fiopenmp', language: ['c', 'cpp'])
  elif effective_compiler_id == 'nvidia_hpc'
    # Mode 1 always uses CPU-only multicore OpenMP (no GPU)
    openmp_c_args = ['-O3', '-mp=gpu,multicore', '-g']
    add_project_link_arguments(['-mp=gpu,multicore'], language: ['c', 'cpp'])
  elif effective_compiler_id == 'clang'
    # Clang with libomp (e.g., from LLVM, Homebrew on macOS, or AOMP)
    openmp_c_args = ['-O3', '-fopenmp', '-g']
    add_project_link_arguments('-fopenmp', language: ['c', 'cpp'])
  else
    # GCC and other compilers.
    # -foffload=disable: installing gcc-offload-nvptx makes GCC automatically
    # attempt nvptx offloading for all -fopenmp links; disable it here so that
    # sw_domain_openmp_ext (CPU-only) does not trigger the nvptx mkoffload pass.
    # The GPU extension (sw_domain_gpu_ext) sets its own -foffload=nvptx-none.
    openmp_c_args = ['-O3', '-fopenmp', '-foffload=disable', '-g']
    add_project_link_arguments(['-fopenmp', '-foffload=disable'], language: ['c', 'cpp'])
  endif
  openmp_deps = dependencies + [openmp]
else
  openmp_c_args = []
  openmp_deps = dependencies
endif


#========================================================
# Install the Python scripts
#========================================================
conf = configuration_data()
conf.set('PYTHON', py3.full_path())


message('BinDir: ' + get_option('bindir'))
configure_file(
    input: 'scripts/anuga_pmesh_gui.py',
    output: 'anuga_pmesh_gui',
    configuration: conf,
    install_dir: get_option('bindir'),
    install_mode: 'rwxr-xr-x'
)

configure_file(
    input: 'scripts/anuga_sww_merge.py',
    output: 'anuga_sww_merge',
    configuration: conf,
    install_dir: get_option('bindir'),
    install_mode: 'rwxr-xr-x'
)

configure_file(
    input: 'scripts/anuga_benchmark_omp.py',
    output: 'anuga_benchmark_omp',
    configuration: conf,
    install_dir: get_option('bindir'),
    install_mode: 'rwxr-xr-x'
)

configure_file(
    input: 'scripts/anuga_toml_run.py',
    output: 'anuga_toml_run',
    configuration: conf,
    install_dir: get_option('bindir'),
    install_mode: 'rwxr-xr-x'
)

configure_file(
    input: 'scripts/anuga_animate_gui.py',
    output: 'anuga_animate_gui',
    configuration: conf,
    install_dir: get_option('bindir'),
    install_mode: 'rwxr-xr-x'
)

configure_file(
    input: 'scripts/anuga_sww_gui.py',
    output: 'anuga_sww_gui',
    configuration: conf,
    install_dir: get_option('bindir'),
    install_mode: 'rwxr-xr-x'
)

configure_file(
    input: 'scripts/anuga_run_isolated_tests.py',
    output: 'anuga_run_isolated_tests',
    configuration: conf,
    install_dir: get_option('bindir'),
    install_mode: 'rwxr-xr-x'
)

# If mpicc is not on meson's PATH (common in conda envs where the env bin
# dir is in the user's PATH but not in the subprocess PATH seen by meson),
# locate it via mpi4py config + shutil.which + sys.prefix/bin fallback, then
# register it so dependency('mpi', ...) in subdirectories finds it correctly.
_mpicc_from_mpi4py = run_command(py3, ['-c', '''
import sys, os, shutil
try:
    if sys.platform == "win32":
        sys.exit(1)  # Windows has no mpicc wrapper; skip override
    import mpi4py
    cfg = mpi4py.get_config()
    mpicc = cfg.get("mpicc", "")
    # Resolve bare name to absolute path
    if mpicc and not os.path.isabs(mpicc):
        mpicc = shutil.which(mpicc) or ""
    # Fallback: look directly in the active Python env's bin/
    if not mpicc or not os.path.isfile(mpicc):
        candidate = os.path.join(sys.prefix, "bin", "mpicc")
        if os.path.isfile(candidate):
            mpicc = candidate
    if mpicc and os.path.isfile(mpicc):
        print(mpicc)
    else:
        sys.exit(1)
except Exception:
    sys.exit(1)
'''], check: false)

if _mpicc_from_mpi4py.returncode() == 0
  _mpicc_path = _mpicc_from_mpi4py.stdout().strip()
  _mpicc_prog = find_program(_mpicc_path, required: false)
  if _mpicc_prog.found()
    meson.override_find_program('mpicc', _mpicc_prog)
    message('mpicc registered from mpi4py: ' + _mpicc_path)
  endif
endif

# Add subdirectories which contains python sources
subdir('anuga')





