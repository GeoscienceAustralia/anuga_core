

inc_dir = include_directories('../utilities', incdir_numpy)

# Compile the Cython-generated C code and additional C code
py3.extension_module('advection_ext',
  sources: ['advection_ext.pyx'],
  include_directories: inc_dir,
  dependencies: dependencies,
  subdir: 'anuga/advection',
  install: true,
)

python_sources = [
  'advection.py', 
  '__init__.py'
  ]

py3.install_sources(
  python_sources,
  subdir: 'anuga/advection'
)

subdir('tests')