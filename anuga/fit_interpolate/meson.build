
inc_dir = include_directories('../utilities')

util_srcs = ['../utilities/quad_tree.c',
                 '../utilities/sparse_dok.c',
                 '../utilities/sparse_csr.c']

# Compile the Cython-generated C code and additional C code
py3.extension_module('fitsmooth_ext',
  sources: ['fitsmooth_ext.pyx'] + util_srcs,
  include_directories: inc_dir,
  dependencies: dependencies,
  subdir: 'anuga/fit_interpolate',
  install: true,
)

python_sources = [
  'benchmark_least_squares.py',
  'fit.py',
  'general_fit_interpolate.py',
  '__init__.py',
  'interpolate2d.py',
  'interpolate.py',
  'profile_long_benchmark.py',
  'run_long_benchmark.py',
  'ticket178_benchmark.py',
]

py3.install_sources(
  python_sources,
  subdir: 'anuga/fit_interpolate'
)

subdir('tests')