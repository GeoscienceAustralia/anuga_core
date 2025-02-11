py3 = import('python').find_installation(pure: false)
dep_py = py3.dependency()

# Find NumPy
dep_numpy = dependency('numpy')

dependencies = [dep_py, dep_numpy]

inc_dir = include_directories('../utilities')

# Compile the Cython-generated C code and additional C code
# py3.extension_module('advection_ext',
#   sources: ['advection_ext.pyx'],
#   include_directories: inc_dir,
#   dependencies: dependencies,
#   subdir: 'anuga/advection',
#   install: true,
# )

python_sources = [
  'AppShell.py',
  'cursornames.py',
  'exesetup.py',
  'graphical_mesh_generator.py',
  '__init__.py',
  'mesh_interface.py',
  'mesh.py',
  'mesh_quadtree.py',
  'ProgressBar.py',
  'setup.py',
  'timing.py',
  'toolbarbutton.py',
  'ungen_example.py',
  'visualmesh.py',
]

py3.install_sources(
  python_sources,
  subdir: 'anuga/pmesh'
)

subdir('tests')