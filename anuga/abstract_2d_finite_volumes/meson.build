
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

# FIXME SR: Need to setup dependencies for openacc code
openacc_deps = dependencies

# Compile the Cython-generated C code and additional C code
py3.extension_module('mesh_factory_ext',
  sources: ['mesh_factory_ext.pyx'],
  include_directories: inc_dir,
  dependencies: dependencies,
  subdir: 'anuga/abstract_2d_finite_volumes',
  install: true,
)

py3.extension_module('neighbour_mesh_ext',
  sources: ['neighbour_mesh_ext.pyx'],
  include_directories: inc_dir,
  dependencies: dependencies,
  subdir: 'anuga/abstract_2d_finite_volumes',
  install: true,
)

py3.extension_module('neighbour_table_ext',
  sources: ['neighbour_table_ext.pyx'],
  override_options : ['cython_language=cpp'],
  include_directories: inc_dir,
  dependencies: dependencies,
  subdir: 'anuga/abstract_2d_finite_volumes',
  install: true,
)

py3.extension_module('pmesh2domain_ext',
  sources: ['pmesh2domain_ext.pyx'],
  include_directories: inc_dir,
  dependencies: dependencies,
  subdir: 'anuga/abstract_2d_finite_volumes',
  install: true,
)

py3.extension_module('quantity_ext',
  sources: ['quantity_ext.pyx'],
  include_directories: inc_dir,
  dependencies: dependencies,
  subdir: 'anuga/abstract_2d_finite_volumes',
  install: true,
)


py3.extension_module('quantity_openmp_ext',
  sources: ['quantity_openmp_ext.pyx'],
  include_directories: inc_dir,
  dependencies: openmp_deps,
  subdir: 'anuga/abstract_2d_finite_volumes',
  install: true,
)

py3.extension_module('quantity_openacc_ext',
  sources: ['quantity_openacc_ext.pyx'],
  include_directories: inc_dir,
  dependencies: openacc_deps,
  subdir: 'anuga/abstract_2d_finite_volumes',
  install: true,
)

py3.extension_module('quantity_cuda_ext',
  sources: ['quantity_cuda_ext.pyx'],
  include_directories: inc_dir,
  dependencies: dependencies,
  subdir: 'anuga/abstract_2d_finite_volumes',
  install: true,
)



python_sources = [
  'ermapper_grids.py',
  'file_function.py',
  'gauge.py',
  'general_mesh.py',
  'generic_boundary_conditions.py',
  'generic_domain.py',
  '__init__.py',
  'mesh_factory.py',
  'neighbour_mesh.py',
  'old_setup.py',
  'pmesh2domain.py',
  'quantity.py',
  'region.py',
  'tag_region.py',
  'util.py',
]


py3.install_sources(
  python_sources,
  subdir: 'anuga/abstract_2d_finite_volumes'
)

subdir('tests')
