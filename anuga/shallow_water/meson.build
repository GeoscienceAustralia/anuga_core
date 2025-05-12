
inc_dir = include_directories('../utilities', incdir_numpy)

openmp = dependency('openmp', required: false)
if openmp.found()
  if host_machine.system() == 'windows'
    # On Windows, the mingw compiler does not support OpenMP ATOMIC operations
    openmp_deps = dependencies
  else
    openmp_deps = dependencies + [openmp]
  endif

else
  openmp_deps = dependencies
endif

# FIXME SR: Need to setup depenencies for openacc code


# py3.extension_module('boundaries',
#   sources: ['boundaries.pyx'],
#   include_directories: inc_dir,
#   dependencies: dependencies,
#   subdir: 'anuga/shallow_water',
#   install: true,
# )

# Compile the Cython-generated C code and additional C code
py3.extension_module('sw_domain_orig_ext',
  sources: ['sw_domain_orig_ext.pyx'],
  include_directories: inc_dir,
  dependencies: dependencies,
  subdir: 'anuga/shallow_water',
  install: true,
)

py3.extension_module('sw_domain_simd_ext',
  sources: ['sw_domain_simd_ext.pyx'],
  include_directories: inc_dir,
  dependencies: dependencies,
  subdir: 'anuga/shallow_water',
  install: true,
)

py3.extension_module('sw_domain_openmp_ext',
  sources: ['sw_domain_openmp_ext.pyx'],
  c_args : ['-O3', '-march=native'],
  include_directories: inc_dir,
  dependencies: openmp_deps,
  subdir: 'anuga/shallow_water',
  install: true,
)


py3.extension_module('sw_domain_openacc_ext',
  sources: ['sw_domain_openacc_ext.pyx'],
  include_directories: inc_dir,
  dependencies: dependencies,
  subdir: 'anuga/shallow_water',
  install: true,
) 



python_sources = [
'boundaries.py',
'checkpoint.py',
'forcing.py',
'friction.py',
'__init__.py',
'most2nc.py',
'shallow_water_domain.py',
'sww_interrogate.py',
'tsh2sww.py',
]


py3.install_sources(
  python_sources,
  subdir: 'anuga/shallow_water'
)

subdir('tests')
