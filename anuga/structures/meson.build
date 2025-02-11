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
'boyd_box_operator_Amended3.py',
'boyd_box_operator.py',
'boyd_pipe_operator.py',
'__init__.py',
'inlet_enquiry.py',
'inlet_operator.py',
'inlet.py',
'internal_boundary_functions.py',
'internal_boundary_operator.py',
'riverwall.py',
'setup.py',
'structure_operator.py',
'weir_orifice_trapezoid_operator.py',
  ]

py3.install_sources(
  python_sources,
  subdir: 'anuga/structures'
)

subdir('tests')