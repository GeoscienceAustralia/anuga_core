
inc_dir = include_directories('../utilities')


# Compile the Cython-generated C code and additional C code
py3.extension_module('sw_domain_orig_ext',
  sources: ['sw_domain_orig_ext.pyx'],
  include_directories: inc_dir,
  dependencies: dependencies,
  subdir: 'anuga/shallow_water',
  install: true,
)


python_sources = [
'boundaries.py',
'checkpoint.py',
'forcing.py',
'__init__.py',
'most2nc.py',
'setup.py',
'shallow_water_domain.py',
'sww_interrogate.py',
'tsh2sww.py',
]


py3.install_sources(
  python_sources,
  subdir: 'anuga/shallow_water'
)

subdir('tests')
