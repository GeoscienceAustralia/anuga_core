
inc_dir = include_directories('../utilities', incdir_numpy)

py3.extension_module('mesh_relax_ext',
  sources: ['mesh_relax_ext.pyx'],
  include_directories: inc_dir,
  dependencies: dependencies,
  subdir: 'anuga/pmesh',
  install: true,
)

python_sources = [
  'AppShell.py',
  'cursornames.py',
  'exesetup.py',
  '__init__.py',
  'mesh_interface.py',
  'mesh.py',
  'mesh_quadtree.py',
  'mesh_relax.py',
  'ProgressBar.py',
  'size_field.py',
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
subdir('icons')
