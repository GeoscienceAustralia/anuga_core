project('anuga', 'c', 'cpp', 'cython', version: '3.2.dev0')

# Import python
py3 = import('python').find_installation(pure: false)

message('Python found at ' + py3.get_variable('prefix'))

dep_py = py3.dependency()

# Find NumPy
dep_numpy = dependency('numpy', method : 'config-tool', required : true)

# important to put the numpy dependency first to avoid picking up the system numpy
dependencies = [dep_numpy, dep_py]

# Add subdirectories which contains python sources
subdir('anuga')





