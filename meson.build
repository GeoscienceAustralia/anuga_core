project('anuga', 'c', 'cpp', 'cython', version: '3.2.dev0', default_options: ['cpp_std=c++17'])

add_project_arguments('-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION', language: 'cpp')
add_project_arguments('-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION', language: 'c')
# if windows
#add_project_arguments('-DMS_WIN64', language: 'cpp')
#add_project_arguments('-DMS_WIN64', language: 'c')

# Import python
py3 = import('python').find_installation(pure: false)

message('Python found at ' + py3.get_variable('prefix'))

dep_py = py3.dependency()

# Find NumPy
dep_numpy = dependency('numpy')


# incdir_numpy = run_command(py3,
#   [
#     '-c',
#     'import os; os.chdir(".."); import numpy; print(numpy.get_include())'
#   ],
#   check: true
# ).stdout().strip()

# #inc_np = include_directories(incdir_numpy)

# #incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
# #inc_f2py = include_directories(incdir_f2py)
# #fortranobject_c = incdir_f2py / 'fortranobject.c'

# npymath_path = incdir_numpy / '..' / 'lib'
# npymath_lib = meson.get_compiler('c').find_library('npymath', dirs: npymath_path)

# my_inc = include_directories(incdir_numpy)
# my_lib = shared_library(npymath_path / 'libnpymath.a')
# my_dep = declare_dependency(link_with : my_lib,
#   include_directories : my_inc)

# important to put the numpy dependency first to avoid picking up the system numpy
dependencies = [dep_numpy, dep_py]
#dependencies = [dep_py]

# Add subdirectories which contains python sources
subdir('anuga')





