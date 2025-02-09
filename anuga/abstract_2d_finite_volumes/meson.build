# meson.build inside anuga/abstract_2d_finite_volumes

# project('anuga', 'cython', version: '1.0')

# # Import python and cython modules
# py3 = import('python').find_installation()
# dep_py = py3.dependency()


# #cython = import('cython')

# # Find NumPy
# dep_numpy = dependency('numpy')

# dependencies = [dep_py, dep_numpy]

# Cython to C conversion
#cython_files = cython.generate(pyx_sources: files('*.pyx'))

inc_dir = include_directories('../utilities')

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
  python_sources
)

subdir('tests')
