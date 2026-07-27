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

# cc = meson.get_compiler('c')
# if openmp_dep.found()
#   if cc.get_id() in ['intel', 'intel-cl', 'icx']
#     openmp_flag = '-qopenmp'
#   else
#     openmp_flag = '-fopenmp'
#   endif
# endif

openmp = dependency('openmp', required: false)
if openmp.found()
  # On Windows, the mingw compiler does not support OpenMP ATOMIC operations
  # so using gcc_win-64 gxx_win-64 provided by conda-forge
  # On linux, OpenMP is supported by gcc and g++
  # On macOS, OpenMP is not supported by the Apple clang compiler so using 
  # clang and llvm-openmp provided by conda-forge
 
  #if host_machine.system() == 'windows'
  
  # Do NOT use -march=native here: it bakes the build host's ISA into the
  # compiled extensions, so a wheel built on e.g. a Sapphire Rapids CI runner
  # gets AVX-512 and dies with an illegal instruction at `import anuga` on any
  # CPU without it (this broke the 3.3.8 PyPI wheels, e.g. on Colab). Leave the
  # toolchain's portable default in place (conda-forge sets -march=nocona).
  openmp_c_args = ['-O3', '-fopenmp', '-g']

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
    input: 'scripts/anuga_run_toml.py',
    output: 'anuga_run_toml',
    configuration: conf,
    install_dir: get_option('bindir'),
    install_mode: 'rwxr-xr-x'
)

# Add subdirectories which contains python sources
subdir('anuga')





