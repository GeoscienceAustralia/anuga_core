
inc_dir = include_directories('../utilities', incdir_numpy)

# Compile the Cython-generated C code and additional C code
py3.extension_module('mannings_operator_ext',
  sources: ['mannings_operator_ext.pyx'],
  include_directories: inc_dir,
  dependencies: dependencies,
  subdir: 'anuga/operators',
  install: true,
)

py3.extension_module('kinematic_viscosity_operator_ext',
  sources: ['kinematic_viscosity_operator_ext.pyx'],
  include_directories: inc_dir,
  dependencies: dependencies,
  subdir: 'anuga/operators',
  install: true,
)

python_sources = [
  'base_operator.py',
  'boundary_flux_integral_operator.py',
  'change_friction_operator.py',
  'collect_max_quantities_operator.py',
  'collect_max_stage_operator.py',
  'elliptic_operator.py',
  'erosion_operators.py',
  '__init__.py',
  'kinematic_viscosity_operator.py',
  'mannings_operator.py',
  'rate_operators.py',
  'sanddune_erosion_operator.py',
  'set_elevation_operator.py',
  'set_elevation.py',
  'set_friction_operators.py',
  'set_quantity_operator.py',
  'set_quantity.py',
  'set_stage_operator.py',
  'set_stage.py',
  'set_w_uh_vh_operator.py',
]

py3.install_sources(
  python_sources,
  subdir: 'anuga/operators'
)

subdir('tests')