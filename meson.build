project('anuga', 'c', 'cpp', 'cython', version: '3.2.dev0')

# Import python
py3 = import('python').find_installation(pure: false, modules: ['numpy'])

message('Python found at ' + py3.get_variable('prefix'))

dep_py = py3.dependency()

# Find NumPy
#dep_numpy = dependency('numpy', required : true)


incdir_numpy = run_command(py3,
  [
    '-c',
    'import os; os.chdir(".."); import numpy; print(numpy.get_include())'
  ],
  check: true
).stdout().strip()

inc_np = include_directories(incdir_numpy)

#incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
#inc_f2py = include_directories(incdir_f2py)
#fortranobject_c = incdir_f2py / 'fortranobject.c'

npymath_path = incdir_numpy / '..' / 'lib'
npymath_lib = meson.get_compiler('c').find_library('npymath', dirs: npymath_path)

my_inc = include_directories(inc_np)
my_lib = static_library(npymath_lib)
my_dep = declare_dependency(link_with : my_lib,
  include_directories : my_inc)

# important to put the numpy dependency first to avoid picking up the system numpy
dependencies = [my_dep, dep_py]

# Add subdirectories which contains python sources
subdir('anuga')





