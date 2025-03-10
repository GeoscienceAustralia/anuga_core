
inc_dir = include_directories('../utilities', incdir_numpy)

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
