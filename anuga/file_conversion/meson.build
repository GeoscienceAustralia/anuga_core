
inc_dir = include_directories('../utilities', incdir_numpy)

# Compile the Cython-generated C code and additional C code
py3.extension_module('calc_grid_values_ext',
  sources: ['calc_grid_values_ext.pyx'],
  include_directories: inc_dir,
  dependencies: dependencies,
  subdir: 'anuga/file_conversion',
  install: true,
)

python_sources = [
  'asc2dem.py',
  'csv2sts.py',
  'dem2array.py',
  'dem2dem.py',
  'dem2pts.py',
  'esri2sww.py',
  'ferret2sww.py',
  'file_conversion.py',
  'grd2array.py',
  '__init__.py',
  'llasc2pts.py',
  'sdf2pts.py',
  'sts2sww_mesh.py',
  'sww2array.py',
  'sww2dem_new.py',
  'sww2dem.py',
  'sww2pts.py',
  'tif2array.py',
  'tif2point_values.py',
  'urs2nc.py',
  'urs2sts.py',
  'urs2sww.py',
  'urs2txt.py',
  'xya2pts.py',
]

py3.install_sources(
  python_sources,
  subdir: 'anuga/file_conversion'
)

subdir('tests')