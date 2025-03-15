project('anuga', 'c', 'cpp', 'cython', version: '3.2.dev0', default_options: ['cpp_std=c++17'])

add_project_arguments('-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION', language: 'cpp')
add_project_arguments('-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION', language: 'c')

py3 = import('python').find_installation(pure: false)

dep_py = py3.dependency()

# Find NumPy
dep_numpy = dependency('numpy')

incdir_numpy = run_command(py3,
  [
    '-c',
    'import os; os.chdir(".."); import numpy; print(numpy.get_include())'
  ],
  check: true
).stdout().strip()

message('NumPy include directory: ' + incdir_numpy)

# important to put the numpy dependency first to avoid picking up the system numpy
dependencies = [dep_py, dep_numpy]

conf = configuration_data()
conf.set('PYTHON', py3.path())

configure_file(
    input: 'anuga/pmesh/graphical_mesh_generator.py',
    output: 'graphical_mesh_generator',
    configuration: conf,
    install_dir: get_option('bindir')
)

# Add subdirectories which contains python sources
subdir('anuga')





