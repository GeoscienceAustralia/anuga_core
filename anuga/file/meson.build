py3 = import('python').find_installation(pure: false)
dep_py = py3.dependency()

# Find NumPy
dep_numpy = dependency('numpy')

dependencies = [dep_py, dep_numpy]

inc_dir = include_directories('../utilities')

# Compile the Cython-generated C code and additional C code
py3.extension_module('urs_ext',
  sources: ['urs_ext.pyx'],
  include_directories: inc_dir,
  dependencies: dependencies,
  subdir: 'anuga/file',
  install: true,
)

python_sources = [
  'csv_file.py',
  '__init__.py',
  'mux.py',
  'netcdf.py',
  'pts.py',
  'sts.py',
  'sww.py',
  'ungenerate.py',
  'urs.py',
]

py3.install_sources(
  python_sources,
  subdir: 'anuga/file'
)

subdir('tests')